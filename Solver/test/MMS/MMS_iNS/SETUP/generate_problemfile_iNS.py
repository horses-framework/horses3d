#!/usr/bin/env python3
"""
generate_problemfile_iNS.py
===========================
Generates ProblemFile.f90 for the HORSES3D Incompressible Navier-Stokes (iNS) MMS convergence study.

USAGE
-----
1. Edit SECTION 1 — Physical parameters (must match control file).
   - RHO0: reference density (should match control)
2. Edit SECTION 2 — Manufactured solution (must be periodic on [-1,1]^3 for periodic BC).
   - Density can be constant (ρ = RHO0) or variable.
   - Velocity and pressure fields must be defined.
3. Edit SECTION 3 — Output settings.
4. Run:  python3 generate_problemfile_iNS.py
5. Compile the generated ProblemFile.f90 with your HORSES3D build.
"""

import sympy as sp
import numpy as np
import sys
import signal
import os

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — Physical parameters
#  These must match the control file.
#  RHO0:   reference density (make sure it matches control)
# ═══════════════════════════════════════════════════════════════════════════════

RHO0 = 1.0        # reference density 

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — Manufactured solution
#  Must be periodic on [-1,1]^3 for periodic BC.
#  Density can be constant or variable; for incompressible flow, constant is typical.
#  Velocity field is chosen divergence-free for consistency, but source terms will
#  enforce the equations if not.
#
#  IMPORTANT: All fields must be SymPy expressions. If you want a constant,
#             use sp.Rational or sp.Number, e.g., rho = sp.Rational(1,1)
# ═══════════════════════════════════════════════════════════════════════════════

t, x, y, z = sp.symbols('t x y z')
π = sp.pi

# Density (constant by default; can be made variable)
rho = RHO0   # constant example; replace with your own expression
# Example variable density:
# rho = RHO0 * (1 + 0.1*sp.sin(π*x)*sp.sin(π*y)*sp.sin(π*z)*sp.cos(t))

# Velocity components
u = sp.sin(π*x)*sp.cos(π*y)*sp.cos(π*z)*sp.sin(t)
v = sp.cos(π*x)*sp.sin(π*y)*sp.cos(π*z)*sp.sin(t)
w = -2*sp.cos(π*x)*sp.cos(π*y)*sp.sin(π*z)*sp.sin(t)

# Pressure (must be positive)
p = 1 + 0.1*sp.sin(π*x)*sp.sin(π*y)*sp.sin(π*z)*sp.cos(t)

# Optional: replace the above with your own expressions.
# Example non‑periodic test:
# u = -sp.cos(x)*sp.sin(y)*sp.sin(z)*sp.exp(-t)
# v =  sp.sin(x)*sp.cos(y)*sp.sin(z)*sp.exp(-t)
# w = -sp.sin(x)*sp.sin(y)*sp.cos(z)*sp.exp(-t)
# p = (sp.cos(x)+sp.cos(y)+sp.cos(z))*sp.exp(-2*t)
# rho = RHO0

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3 — Output settings
# ═══════════════════════════════════════════════════════════════════════════════

OUTPUT_FILE = "ProblemFile.f90"

# Quadrature node type — must match "Discretization nodes" in control file
# "Gauss-Lobatto"  or  "Gauss"
NODE_TYPE = "Gauss-Lobatto"

# Simplification strategy for SymPy source terms:
# "auto"     — try simplify(), fall back to trigsimp after TIMEOUT seconds
# "simplify" — always use simplify() (slow, works always)
# "trigsimp" — always use trigsimp() (fast; use on Windows or for pure trig)
SIMPLIFY_STRATEGY = "auto"
SIMPLIFY_TIMEOUT  = 60

# ═══════════════════════════════════════════════════════════════════════════════
#  END OF USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# ── Convert numeric fields to SymPy expressions ────────────────────────────────
def to_sympy(expr):
    if isinstance(expr, (int, float)):
        return sp.Number(expr)
    return expr

rho = to_sympy(rho)
u   = to_sympy(u)
v   = to_sympy(v)
w   = to_sympy(w)
p   = to_sympy(p)

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

print("Computing iNS source terms...", file=sys.stderr)

# Symbols for physical parameters
rho_0_sym = sp.Symbol('rho_0')       # reference density
c0_sq_sym = sp.Symbol('c0_sq')      # artificial speed of sound squared (thermodynamics_ % rho0c02)
mu_sym    = sp.Symbol('mu')          # dynamic viscosity

# Manufactured fields

# Velocity gradients
ux, uy, uz = sp.diff(u,x), sp.diff(u,y), sp.diff(u,z)
vx, vy, vz = sp.diff(v,x), sp.diff(v,y), sp.diff(v,z)
wx, wy, wz = sp.diff(w,x), sp.diff(w,y), sp.diff(w,z)

# Divergence of velocity
divU = ux + vy + wz

# Viscous stress tensor components (incompressible form)
tau_xx = mu_sym * (2*ux)            # no bulk viscosity term
tau_yy = mu_sym * (2*vy)
tau_zz = mu_sym * (2*wz)
tau_xy = mu_sym * (uy + vx)
tau_xz = mu_sym * (uz + wx)
tau_yz = mu_sym * (vz + wy)

# Divergence of viscous stress (for momentum equations)
div_tau_x = div(tau_xx, tau_xy, tau_xz)
div_tau_y = div(tau_xy, tau_yy, tau_yz)
div_tau_z = div(tau_xz, tau_yz, tau_zz)

# RHS expressions for continuity, momentum, pressure
S_rho   = sp.diff(rho, t) + div(rho*u, rho*v, rho*w)
S_rhou  = sp.diff(rho*u, t) + div(rho*u*u + p, rho*u*v, rho*u*w) - div_tau_x
S_rhov  = sp.diff(rho*v, t) + div(rho*v*u, rho*v*v + p, rho*v*w) - div_tau_y
S_rhow  = sp.diff(rho*w, t) + div(rho*w*u, rho*w*v, rho*w*w + p) - div_tau_z
S_p     = sp.diff(p, t) + c0_sq_sym * divU

rhs_exprs = [S_rho, S_rhou, S_rhov, S_rhow, S_p]
var_names = ['INSRHO', 'INSRHOU', 'INSRHOV', 'INSRHOW', 'INSP']

rhs_simplified = []
for i, expr in enumerate(rhs_exprs):
    print(f"  S({var_names[i]}) [{i+1}/{len(rhs_exprs)}]...", end=' ', flush=True, file=sys.stderr)
    rhs_simplified.append(simp(expr))
    print("done", file=sys.stderr)

# ── Fortran code emission ──────────────────────────────────────────────────────

PI_sym = sp.Symbol('PI')
x1     = sp.Symbol('x(1)')
x2     = sp.Symbol('x(2)')
x3     = sp.Symbol('x(3)')

def to_f90(expr, lhs, indent='            '):
    # Replace SymPy symbols with Fortran identifiers
    s = expr.subs(sp.pi, PI_sym)
    s = s.subs(x, x1).subs(y, x2).subs(z, x3)
    s = s.subs(rho_0_sym, sp.Symbol('rho_0'))
    s = s.subs(c0_sq_sym, sp.Symbol('c0_sq'))
    s = s.subs(mu_sym, sp.Symbol('mu'))
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
    # Substitute the reference density symbol with the Fortran parameter RHO0
    rho_expr = rho.subs(rho_0_sym, sp.Symbol('RHO0'))
    u_expr   = u.subs(sp.pi, PI_sym)
    v_expr   = v.subs(sp.pi, PI_sym)
    w_expr   = w.subs(sp.pi, PI_sym)
    p_expr   = p.subs(sp.pi, PI_sym)

    # Conserved variables: rho, rhou, rhov, rhow, p
    res_rho  = sp.fcode(rho_expr, assign_to='res(INSRHO)',  source_format='free', standard=95)
    res_rhou = sp.fcode(rho_expr*u_expr, assign_to='res(INSRHOU)', source_format='free', standard=95)
    res_rhov = sp.fcode(rho_expr*v_expr, assign_to='res(INSRHOV)', source_format='free', standard=95)
    res_rhow = sp.fcode(rho_expr*w_expr, assign_to='res(INSRHOW)', source_format='free', standard=95)
    res_p    = sp.fcode(p_expr, assign_to='res(INSP)', source_format='free', standard=95)

    # Indent and combine
    def indent_block(code, ind='      '):
        return '\n'.join(ind + l for l in code.split('\n'))

    return f"""\
   subroutine evalManufacturedSolution(x, y, z, t, res)
!
!     Evaluates the manufactured solution at (x,y,z,t).
!     Conserved variables: rho, rhou, rhov, rhow, p
!
!     RHO0 : reference density 
!
      use SMConstants
      use PhysicsStorage_iNS
      implicit none
      real(kind=RP), intent(in)  :: x, y, z, t
      real(kind=RP), intent(out) :: res(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP), parameter :: RHO0 = {RHO0}_RP

{indent_block(res_rho)}
{indent_block(res_rhou)}
{indent_block(res_rhov)}
{indent_block(res_rhow)}
{indent_block(res_p)}

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
!      ProblemFile.f90  —  HORSES3D Incompressible Navier-Stokes (iNS) MMS Study
!      Generated by generate_problemfile_iNS.py
!
!      Manufactured solution (periodic on [-1,1]^3):
!        rho = defined by user
!        u = sin(pi*x)*cos(pi*y)*cos(pi*z)*sin(t)
!        v = cos(pi*x)*sin(pi*y)*cos(pi*z)*sin(t)
!        w = -2*cos(pi*x)*cos(pi*y)*sin(pi*z)*sin(t)
!        p = 1 + 0.1*sin(pi*x)*sin(pi*y)*sin(pi*z)*cos(t)
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

#ifdef INCNS
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
            real(kind=RP) :: t, rho_0, mu, c0_sq

#ifdef INCNS
            t     = time
            ! rho_0: reference density 
            rho_0 = {RHO0}_RP
            ! mu: dynamic viscosity (isotropic, take first component)
            mu    = dimensionless_ % mu(1)
            ! c0_sq: artificial speed of sound squared (thermodynamics_ % rho0c02)
            c0_sq = thermodynamics_ % rho0c02

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

#ifdef INCNS
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
            print *, "INCNS MMS L2 error: ", error_mesh
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
print(f"Next step: compile with 'make' in your HORSES3D setup directory.", file=sys.stderr)
