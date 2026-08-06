"""
generate_problemfile_MU.py
==========================
Generates ProblemFile.f90 for the HORSES3D Multiphase (MU) MMS convergence study.

USAGE
-----
1. Edit SECTION 1 — Physical parameters (must match control file).
2. Edit SECTION 2 — Manufactured solution (must be periodic on [-1,1]^3).
3. Edit SECTION 3 — Output settings.
4. Run:  python3 generate_problemfile_MU.py
5. Compile the generated ProblemFile.f90, it is not automatic

The generated file is self-contained and can be reused across runs as long
as the physical parameters and manufactured solution do not change.

To edit the problem file that will result go to lines 312 onwards

Equations can be changed at line 188 (to add gravity for example)
"""

import sympy as sp
import sys
import signal
import os

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — Physical parameters
#  These must match control file.
#  RHO1 = fluid 1 density / reference density  
#  RHO2 = fluid 2 density / reference density  (= fluid_2_density / fluid_1_density)
# ═══════════════════════════════════════════════════════════════════════════════

RHO1  = 1.0    # dimensionless_ % rho(1)  — adjust to match control file
RHO2  = 2.0    # dimensionless_ % rho(2)  — adjust to match control file

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — Manufactured solution
#  Must be periodic on [-1,1]^3:  sin(pi*(±1)) = 0, cos(pi*(±1)) = -1
#  Velocity field chosen to be divergence-free: du/dx + dv/dy + dw/dz = 0
# ═══════════════════════════════════════════════════════════════════════════════

t, x, y, z = sp.symbols('t x y z')
π = sp.pi

c =  sp.Rational(1,2)*(1 + sp.cos(π*x)*sp.cos(π*y)*sp.cos(π*z)*sp.sin(t))
u =  sp.sin(π*x)*sp.cos(π*y)*sp.cos(π*z)*sp.sin(t)
v =  sp.cos(π*x)*sp.sin(π*y)*sp.cos(π*z)*sp.sin(t)
w = -2*sp.cos(π*x)*sp.cos(π*y)*sp.sin(π*z)*sp.sin(t)
p =  sp.sin(π*x)*sp.sin(π*y)*sp.sin(π*z)*sp.cos(t)

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3 — Output settings
# ═══════════════════════════════════════════════════════════════════════════════

OUTPUT_FILE = "ProblemFile.f90"

# Quadrature node type — must match "Discretization nodes" in control file
# "Gauss-Lobatto"  or  "Gauss"
NODE_TYPE = "Gauss-Lobatto"

# Simplification strategy for SymPy source terms:
# "auto"     — try simplify(), fall back to trigsimp() after TIMEOUT seconds
# "simplify" — always use simplify() (slow, works always)
# "trigsimp" — always use trigsimp() (fast; use on Windows or for pure trig)
SIMPLIFY_STRATEGY = "auto"
SIMPLIFY_TIMEOUT  = 60

# ═══════════════════════════════════════════════════════════════════════════════
#  END OF USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# ── Simplification ─────────────────────────────────────────────────────────────

_SIGALRM_AVAILABLE = hasattr(signal, 'SIGALRM')

def _simp_timeout(expr):
    """Run sp.simplify() with a wall-clock timeout.

    On POSIX systems (Linux, macOS) SIGALRM is used for a hard timeout.
    On Windows (and any platform that lacks SIGALRM) we fall back directly
    to trigsimp so the script still runs correctly.
    """
    if not _SIGALRM_AVAILABLE:
        # Windows / non-POSIX: SIGALRM is unavailable — use trigsimp directly.
        print("  SIGALRM not available (Windows?) — using trigsimp directly",
              file=sys.stderr)
        return sp.trigsimp(sp.expand(expr))

    result = [None]
    def _handler(signum, frame): raise TimeoutError()
    try:
        signal.signal(signal.SIGALRM, _handler)
        signal.alarm(SIMPLIFY_TIMEOUT)
        result[0] = sp.simplify(expr)
        signal.alarm(0)
    except TimeoutError:
        print("  simplify timed out — falling back to trigsimp", file=sys.stderr)
        result[0] = sp.trigsimp(sp.expand(expr))
    except Exception:
        signal.alarm(0)
        raise
    return result[0]

def simp(expr):
    if SIMPLIFY_STRATEGY == "simplify":   return sp.simplify(expr)
    elif SIMPLIFY_STRATEGY == "trigsimp": return sp.trigsimp(sp.expand(expr))
    else:                                 return _simp_timeout(expr)

# ── Differential operators ─────────────────────────────────────────────────────

def lap(f):
    return sp.diff(f,x,2) + sp.diff(f,y,2) + sp.diff(f,z,2)

def div(fx, fy, fz):
    return sp.diff(fx,x) + sp.diff(fy,y) + sp.diff(fz,z)

# ── Gauss-Lobatto weights ──────────────────────────────────────────────────────

import numpy as np

def gl_weights(N):
    if N == 0: return np.array([2.0])
    if N == 1: return np.array([1.0, 1.0])
    nodes = np.zeros(N + 1)
    nodes[0] = -1.0;  nodes[-1] = 1.0
    for i in range(1, N):
        nodes[i] = -np.cos(np.pi * i / N)
    for _ in range(100):
        for i in range(1, N):
            xi = nodes[i]
            p0, p1, dp0, dp1 = 1.0, xi, 0.0, 1.0
            for k in range(1, N):
                p2  = ((2*k+1)*xi*p1 - k*p0) / (k+1)
                dp2 = ((2*k+1)*(p1 + xi*dp1) - k*dp0) / (k+1)
                p0, p1, dp0, dp1 = p1, p2, dp1, dp2
            denom = 1.0 - xi**2
            ddp1 = (2*xi*dp1 - N*(N+1)*p1) / denom if abs(denom) > 1e-14 else 0.0
            if abs(ddp1) > 1e-14: nodes[i] -= dp1 / ddp1
    weights = np.zeros(N + 1)
    for i in range(N + 1):
        xi = nodes[i]
        p0, p1 = 1.0, xi
        for k in range(1, N):
            p2 = ((2*k+1)*xi*p1 - k*p0) / (k+1)
            p0, p1 = p1, p2
        weights[i] = 2.0 / (N*(N+1) * p1**2)
    weights[0] = weights[-1] = 2.0 / (N*(N+1))
    return weights

def gauss_weights(N):
    n = N + 1
    weights = np.zeros(n)
    nodes   = np.zeros(n)
    for i in range(n // 2 + 1):
        xi = np.cos(np.pi * (i + 0.75) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, xi
            for k in range(1, n):
                p2 = ((2*k+1)*xi*p1 - k*p0) / (k+1)
                p0, p1 = p1, p2
            dp = n * (p0 - xi*p1) / (1.0 - xi**2)
            dx = p1 / dp
            xi -= dx
            if abs(dx) < 1e-15: break
        nodes[i] = -xi;  nodes[n-1-i] = xi
        w = 2.0 / ((1.0 - xi**2) * dp**2)
        weights[i] = w;  weights[n-1-i] = w
    return weights

MAX_ORDER = 8
wfn       = gl_weights if NODE_TYPE == "Gauss-Lobatto" else gauss_weights
wtag      = "wGL"      if NODE_TYPE == "Gauss-Lobatto" else "wGG"

# ── Source term derivation ─────────────────────────────────────────────────────

print("Computing MULTIPHASE source terms...", file=sys.stderr)

ρ_1, ρ_2 = sp.symbols('rho_1 rho_2')
M0, σ, ϵ = sp.symbols('M0 sigma eps')
η        = sp.symbols('eta')
cs       = sp.symbols('cs')

ρ_mix = (1 - c)*ρ_1 + c*ρ_2
sqrtρ = sp.sqrt(ρ_mix)

dfdc = (σ/ϵ)*(12*c**2*(2*c-2) + 24*c*(1-c)**2)
μ_CH = dfdc - sp.Rational(3,2)*σ*ϵ*lap(c)

ux,uy,uz = sp.diff(u,x), sp.diff(u,y), sp.diff(u,z)
vx,vy,vz = sp.diff(v,x), sp.diff(v,y), sp.diff(v,z)
wx,wy,wz = sp.diff(w,x), sp.diff(w,y), sp.diff(w,z)

τx = div(η*(2*ux), η*(uy+vx), η*(uz+wx))
τy = div(η*(vx+uy), η*(2*vy), η*(vz+wy))
τz = div(η*(wx+uz), η*(wy+vz), η*(2*wz))

adv_u = u*ux + v*uy + w*uz
adv_v = u*vx + v*vy + w*vz
adv_w = u*wx + v*wy + w*wz

rhs_exprs = [
    sp.diff(c,t) + div(c*u, c*v, c*w) - M0*lap(μ_CH),
    sqrtρ*sp.diff(sqrtρ*u,t) + sp.Rational(1,2)*ρ_mix*adv_u
        + c*sp.diff(μ_CH,x) + sp.diff(p,x) - τx,
    sqrtρ*sp.diff(sqrtρ*v,t) + sp.Rational(1,2)*ρ_mix*adv_v
        + c*sp.diff(μ_CH,y) + sp.diff(p,y) - τy,
    sqrtρ*sp.diff(sqrtρ*w,t) + sp.Rational(1,2)*ρ_mix*adv_w
        + c*sp.diff(μ_CH,z) + sp.diff(p,z) - τz,
    sp.diff(p,t) + ρ_mix*cs**2*div(u,v,w),
]

var_names = ['IMC', 'IMSQRHOU', 'IMSQRHOV', 'IMSQRHOW', 'IMP']

rhs_simplified = []
for i, expr in enumerate(rhs_exprs):
    print(f"  S({var_names[i]}) [{i+1}/{len(rhs_exprs)}]...",
          end=' ', flush=True, file=sys.stderr)
    rhs_simplified.append(simp(expr))
    print("done", file=sys.stderr)

# ── Fortran code emission ──────────────────────────────────────────────────────

PI_sym = sp.Symbol('PI')
x1     = sp.Symbol('x(1)')
x2     = sp.Symbol('x(2)')
x3     = sp.Symbol('x(3)')

def to_f90(expr, lhs, indent='            '):
    s   = expr.subs(sp.pi, PI_sym).subs(x, x1).subs(y, x2).subs(z, x3)
    raw = sp.fcode(s, assign_to=lhs, source_format='free',
                   standard=95, contract=False)
    return '\n'.join(indent + line for line in raw.split('\n'))

def weight_arrays_f90(indent='            '):
    lines = [f'{indent}! {NODE_TYPE} quadrature weights, orders 1..{MAX_ORDER}']
    for N in range(1, MAX_ORDER+1):
        w    = wfn(N)
        vals = ', &\n'.join(f'{indent}      {wi:.17e}_RP' for wi in w)
        lines.append(f'{indent}real(kind=RP) :: {wtag}_{N}(0:{N})')
        lines.append(f'{indent}data {wtag}_{N} / &')
        lines.append(f'{vals} /')
        lines.append('')
    return '\n'.join(lines)

def weight_select_f90(order_var, pt_var, res_var, indent='                        '):
    lines = [f'{indent}select case({order_var})']
    for N in range(1, MAX_ORDER+1):
        lines.append(f'{indent}case({N}) ; {res_var} = {wtag}_{N}({pt_var})')
    lines.append(f'{indent}case default')
    lines.append(f'{indent}   print *, "MMS: order", {order_var}, '
                 f'"not tabulated (max={MAX_ORDER})" ; error stop')
    lines.append(f'{indent}end select')
    return '\n'.join(lines)

def eval_sol_f90():
    """Emit evalManufacturedSolution using the user's solution fields."""
    c_expr   = c.subs(sp.pi, PI_sym)
    rho_expr = ((1 - c)*sp.Symbol('RHO1') + c*sp.Symbol('RHO2')).subs(sp.pi, PI_sym)

    def expr_f90(expr, lhs):
        s = expr.subs(sp.pi, PI_sym)
        return sp.fcode(s, assign_to=lhs, source_format='free',
                        standard=95, contract=False)

    c_code   = expr_f90(c,   'c')
    rho_code = expr_f90((1-c)*sp.Symbol('RHO1') + c*sp.Symbol('RHO2'), 'rho')

    sqrho = sp.sqrt((1-c)*sp.Symbol('RHO1') + c*sp.Symbol('RHO2'))

    res_imc    = expr_f90(c,        'res(IMC)')
    res_sqrhou = expr_f90(sqrho*u,  'res(IMSQRHOU)')
    res_sqrhov = expr_f90(sqrho*v,  'res(IMSQRHOV)')
    res_sqrhow = expr_f90(sqrho*w,  'res(IMSQRHOW)')
    res_imp    = expr_f90(p,        'res(IMP)')

    def indent_block(code, ind='      '):
        return '\n'.join(ind + l for l in code.split('\n'))

    return f"""\
   subroutine evalManufacturedSolution(x, y, z, t, res)
!
!     Evaluates the manufactured solution at (x,y,z,t).
!     Conserved variables: c, sqrt(rho)*u/v/w, p
!
!     RHO1, RHO2: dimensionless densities — must match dimensionless_ % rho(1/2).
!
      use SMConstants
      use PhysicsStorage_MU
      implicit none
      real(kind=RP), intent(in)  :: x, y, z, t
      real(kind=RP), intent(out) :: res(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP) :: c, rho
      real(kind=RP), parameter :: RHO1 = {RHO1}_RP
      real(kind=RP), parameter :: RHO2 = {RHO2}_RP

{indent_block(c_code)}
{indent_block(rho_code)}
{indent_block(res_imc)}
{indent_block(res_sqrhou)}
{indent_block(res_sqrhov)}
{indent_block(res_sqrhow)}
{indent_block(res_imp)}

   end subroutine evalManufacturedSolution"""

def source_terms_f90():
    lines = []
    for name, expr in zip(var_names, rhs_simplified):
        lines.append(f'            ! --- S({name}) ---')
        lines.append(to_f90(expr, f'S({name})'))
        lines.append('')
    return '\n'.join(lines)

# ── Assemble the full ProblemFile ──────────────────────────────────────────────

print("Assembling ProblemFile...", file=sys.stderr)

pf = f"""\
!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90  —  HORSES3D Multiphase (MU) MMS Convergence Study
!      Generated by generate_problemfile_MU.py
!
!      Manufactured solution (periodic on [-1,1]^3):
!
!
!      Node type : {NODE_TYPE}
!      L2 error is written to mms_l2_error.dat on completion.
!      Format: nelems, P, NDOF, L2_error, t_final
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f

      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#ifdef FLOW
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         use VariableConversion, only: GetGradientValues_f
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(inout)       :: U(NGRAD)
         procedure(GetGradientValues_f)     :: GetGradients
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedGradVars_f

      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)    :: x(NDIM)
         real(kind=RP), intent(in)    :: t
         real(kind=RP), intent(in)    :: nHat(NDIM)
         real(kind=RP), intent(in)    :: Q(NCONS)
         real(kind=RP), intent(in)    :: U_x(NGRAD)
         real(kind=RP), intent(in)    :: U_y(NGRAD)
         real(kind=RP), intent(in)    :: U_z(NGRAD)
         real(kind=RP), intent(inout) :: flux(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedNeumann_f
#endif

      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors)
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
      END SUBROUTINE UserDefinedPeriodicOperation_f

#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)    :: x(NDIM)
         real(kind=RP),             intent(in)    :: Q(NCONS)
         real(kind=RP),             intent(in)    :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t),    intent(in)    :: thermodynamics_
         type(Dimensionless_t),     intent(in)    :: dimensionless_
         type(RefValues_t),         intent(in)    :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),        intent(in)    :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif

      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                 , thermodynamics_ &
                                                 , dimensionless_  &
                                                 , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                                 , multiphase_ &
#endif
                                                 , monitors, &
                                                   elapsedTime, &
                                                   CPUTime   )
         use SMConstants
         USE HexMeshClass
         use FluidData
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)                        :: mesh
         REAL(KIND=RP)                         :: time
         integer                               :: iter
         real(kind=RP)                         :: maxResidual
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         type(Monitor_t),        intent(in)    :: monitors
         real(kind=RP),          intent(in)    :: elapsedTime
         real(kind=RP),          intent(in)    :: CPUTime
      END SUBROUTINE UserDefinedFinalize_f

      SUBROUTINE UserDefinedTermination_f
         implicit none
      END SUBROUTINE UserDefinedTermination_f
   end interface

end module ProblemFileFunctions
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE UserDefinedStartup
            IMPLICIT NONE
         END SUBROUTINE UserDefinedStartup
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE UserDefinedFinalSetup(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            implicit none
            class(hexmesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer       :: eID, i, j, k, Nx, Ny, Nz
            real(kind=RP) :: x(NDIM)
            real(kind=RP) :: t

            t = 0.0_RP

#ifdef MULTIPHASE
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
                  call evalManufacturedSolution(x(1), x(2), x(3), t, &
                          mesh % elements(eID) % storage % Q(:,i,j,k))
               end do ; end do ; end do
               end associate
            end do
#endif

         end subroutine UserDefinedInitialCondition
!
!////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in)  :: monitors
         END SUBROUTINE UserDefinedPeriodicOperation
!
!////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
#endif
)
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),          intent(in)    :: x(NDIM)
            real(kind=RP),          intent(in)    :: Q(NCONS)
            real(kind=RP),          intent(in)    :: time
            real(kind=RP),          intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            real(kind=RP) :: t
            real(kind=RP) :: rho_1, rho_2, eta, cs, M0, sigma, eps

#ifdef MULTIPHASE
            t     = time
            rho_1 = dimensionless_ % rho(1)
            rho_2 = dimensionless_ % rho(2)
            eta   = dimensionless_ % mu(1)
            cs    = sqrt(thermodynamics_ % c02(1))
            M0    = multiphase_ % M0
            sigma = multiphase_ % sigma
            eps   = multiphase_ % eps

{source_terms_f90()}
#endif

         end subroutine UserDefinedSourceTermNS
#endif
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),          intent(in)    :: elapsedTime
            real(kind=RP),          intent(in)    :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
            INTEGER                            :: eID, i, j, k, Nx, Ny, Nz
            INTEGER                            :: fid, nelems, P, NDOF
            real(kind=RP)                      :: x(NDIM)
            real(kind=RP)                      :: u_ex(NCONS), u_h(NCONS)
            real(kind=RP)                      :: wi, wj, wk
            real(kind=RP)                      :: error_elem, error_mesh
            !
            ! {NODE_TYPE} quadrature weights, orders 1..{MAX_ORDER}
            !
{weight_arrays_f90()}
            error_mesh = 0.0_RP

#ifdef MULTIPHASE
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)
               error_elem = 0.0_RP

               DO k = 0, Nz ; DO j = 0, Ny ; DO i = 0, Nx
                  x   = mesh % elements(eID) % geom % x(:,i,j,k)
                  u_h = mesh % elements(eID) % storage % Q(:,i,j,k)
                  call evalManufacturedSolution(x(1), x(2), x(3), time, u_ex)

{weight_select_f90('Nx', 'i', 'wi')}
{weight_select_f90('Ny', 'j', 'wj')}
{weight_select_f90('Nz', 'k', 'wk')}
                  error_elem = error_elem + wi*wj*wk &
                             * mesh % elements(eID) % geom % jacobian(i,j,k) &
                             * norm2(u_ex - u_h)**2
               END DO ; END DO ; END DO
               error_mesh = error_mesh + error_elem
            END DO

            error_mesh = sqrt(error_mesh / SIZE(mesh % elements))

            ! --- Write L2 error to file for automated convergence study ---
            ! Format: nelems, P, NDOF, L2_error, t_final
            nelems = SIZE(mesh % elements)
            P      = mesh % elements(1) % Nxyz(1)
            NDOF   = nelems * (P+1)**3
            print *, "MULTIPHASE MMS L2 error: ", error_mesh
            open(newunit=fid, file="mms_l2_error.dat", status="replace")
            write(fid, '(I10, A, I5, A, I12, A, ES24.16, A, ES14.6)') &
                nelems, ', ', P, ', ', NDOF, ', ', error_mesh, ', ', time
            close(fid)
#endif

         END SUBROUTINE UserDefinedFinalize
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE UserDefinedTermination
         IMPLICIT NONE
      END SUBROUTINE UserDefinedTermination
!
!////////////////////////////////////////////////////////////////////////
!
{eval_sol_f90()}
"""

# ── Write output ───────────────────────────────────────────────────────────────

# Resolve OUTPUT_FILE relative to this script's own directory so the generator
# can be launched from any working directory (e.g. python3 path/to/generate_*.py).
# Absolute paths are honoured as-is.
if os.path.isabs(OUTPUT_FILE):
    output_path = OUTPUT_FILE
else:
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               OUTPUT_FILE)

with open(output_path, 'w') as f:
    f.write(pf)

print(f"Written {output_path} ({len(pf.splitlines())} lines)", file=sys.stderr)
print(f"Next step: compile with make in your SETUP directory.", file=sys.stderr)
