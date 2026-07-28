#!/usr/bin/env python3
"""
generate_problemfile_NS.py
==========================
Generates ProblemFile.f90 for the HORSES3D Navier-Stokes (NS) MMS convergence study.

USAGE
-----
1. Edit SECTION 1 — Physical parameters (must match control file).
   RHO0 is the reference density (should be 1.0 in nondimensional form).
2. Edit SECTION 2 — Manufactured solution (must be periodic on [-1,1]^3 for periodic BC).
   Density can be variable; use rho_0 (the symbol) for the reference density.
3. Edit SECTION 3 — Output settings.
4. Run:  python3 generate_problemfile_NS.py
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
#  RHO0:   reference density 
#  GAMMA:  ratio of specific heats
#  MAKE SURE THEY MATCH SOLVER/CONTROL FILE
# ═══════════════════════════════════════════════════════════════════════════════

RHO0  = 1.0        # reference density 
GAMMA = 1.4        # specific heat ratio

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — Manufactured solution
#  Must be periodic on [-1,1]^3 for periodic BC.
#  Velocity field is chosen divergence-free: du/dx + dv/dy + dw/dz = 0
#  Pressure is offset to remain positive.
#  Density can be constant or variable.
#
#  IMPORTANT: All fields must be SymPy expressions. If you want a constant,
#             use sp.Rational or sp.Number, e.g., rho = sp.Rational(1,1)
#             or simply rho = RHO0 (which is a float – we will convert it).
# ═══════════════════════════════════════════════════════════════════════════════

t, x, y, z = sp.symbols('t x y z')
π = sp.pi

# Unsteady 3D manufactured solution for the rotating-cylinder study.
# Annulus in x-y (r in [0.5,1]), axis z, z in [-1,1]; z-frequency 1/2 (half-wave
# over the extent, resolved by the 4/6/8 z-layers). Genuinely time-dependent so
# the time-accurate scheme must track the evolving physics on the rotating mesh.
# rho, p > 0 everywhere.
rho = 1.0 + 0.15*sp.sin(π*x)*sp.cos(π*y)*sp.sin(π*z/2)*sp.cos(t)
u   = 1.0 + 0.10*sp.sin(π*x)*sp.cos(π*y)*sp.cos(π*z/2)*sp.cos(t)
v   =       0.10*sp.cos(π*x)*sp.sin(π*y)*sp.cos(π*z/2)*sp.sin(t)
w   =       0.05*sp.cos(π*x)*sp.cos(π*y)*sp.sin(π*z/2)*sp.sin(t)
p   = 1.0 + 0.15*sp.cos(π*x)*sp.sin(π*y)*sp.cos(π*z/2)*sp.cos(t)

# Optional: you can replace the above with your own expressions.
# Example non‑periodic test:
# u = -sp.cos(x)*sp.sin(y)*sp.sin(z)*sp.exp(-t)
# v =  sp.sin(x)*sp.cos(y)*sp.sin(z)*sp.exp(-t)
# w = -sp.sin(x)*sp.sin(y)*sp.cos(z)*sp.exp(-t)
# p = (sp.cos(x)+sp.cos(y)+sp.cos(z))*sp.exp(-2*t)
# rho = RHO0 * (1 + 0.1*sp.sin(π*x)*sp.sin(π*y)*sp.sin(π*z)*sp.cos(t))

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3 — Output settings
# ═══════════════════════════════════════════════════════════════════════════════

OUTPUT_FILE = "ProblemFile.f90"

# Quadrature node type — must match "Discretization nodes" in control file
# "Gauss-Lobatto"  or  "Gauss"
NODE_TYPE = "Gauss"

# Simplification strategy for SymPy source terms. Simplification is NOT required for
# correctness — fcode emits mathematically equivalent Fortran either way — it only makes
# the generated expressions more compact. With the Sutherland (sqrt) viscosity terms,
# simplify() is very slow, so "none" is the default.
# "none"     — emit the raw expression, no simplification (fast; recommended)
# "auto"     — try simplify(), fall back to trigsimp after TIMEOUT seconds
# "simplify" — always use simplify() (slow, works always)
# "trigsimp" — always use trigsimp() (fast; use on Windows or for pure trig)
SIMPLIFY_STRATEGY = "none"
SIMPLIFY_TIMEOUT  = 60

# ═══════════════════════════════════════════════════════════════════════════════
#  END OF USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# ── Convert numeric fields to SymPy expressions ────────────────────────────────
# This ensures that .subs works later even if the user provided a float.
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

def _simp_timeout(expr):
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
    if SIMPLIFY_STRATEGY == "none":       return expr
    elif SIMPLIFY_STRATEGY == "simplify": return sp.simplify(expr)
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

MAX_ORDER = 12
wfn       = gl_weights if NODE_TYPE == "Gauss-Lobatto" else gauss_weights
wtag      = "wGL"      if NODE_TYPE == "Gauss-Lobatto" else "wGG"

# ── Source term derivation ─────────────────────────────────────────────────────

print("Computing NAVIERSTOKES source terms...", file=sys.stderr)

# Symbols for physical parameters (to be replaced by local Fortran variables)
rho_0_sym = sp.Symbol('rho_0')   # reference density (should be 1.0)
gamma_sym = sp.Symbol('gamma')
eta_sym   = sp.Symbol('eta')     # reference dynamic viscosity = dimensionless % mu = 1/Re
Pr_sym    = sp.Symbol('Pr')      # Prandtl number
gammaM2_sym  = sp.Symbol('gammaM2')    # gamma * Mach^2  (= dimensionless % gammaM2)
sdivTref_sym = sp.Symbol('sdivTref')   # S_Sutherland / TRef_Sutherland
Tnorm_sym    = sp.Symbol('Tnorm')      # refValues%T / TRef_Sutherland (= 1 by default)

# Manufactured fields
# The user-defined rho may contain rho_0_sym. We'll treat it as the reference density.

# Total energy per unit volume
e_tot = p/(gamma_sym-1) + sp.Rational(1,2)*rho*(u**2 + v**2 + w**2)

# Velocity gradients
ux, uy, uz = sp.diff(u,x), sp.diff(u,y), sp.diff(u,z)
vx, vy, vz = sp.diff(v,x), sp.diff(v,y), sp.diff(v,z)
wx, wy, wz = sp.diff(w,x), sp.diff(w,y), sp.diff(w,z)

divU = ux + vy + wz

# ── Temperature-dependent viscosity: Sutherland's law ──────────────────────────
# horses3d ALWAYS applies Sutherland's law (VariableConversion_NS.f90,
# get_laminar_mu_kappa) — there is no constant-viscosity switch — so a consistent
# manufactured source MUST use the same variable viscosity/conductivity:
#     mu(T)    = eta * suther(T),        eta   = dimensionless % mu = 1/Re
#     kappa(T) = mu(T) * mu_to_kappa,    mu_to_kappa = 1/((gamma-1) Mach^2 Pr)
# with horses3d's nondimensional temperature  T_h = gammaM2 * p/rho  (Temperature(Q))
# and  suther(T) = (1 + S~)/(T~ + S~) * T~^(3/2),   T~ = Tnorm * T_h.
# The energy viscous flux is  kappa(T)*nablaT  with  nablaT = gammaM2*grad(p/rho)
# (ViscousFlux_STATE), so writing Fourier's law with the manufactured T = p/rho the
# heat-flux coefficient becomes  kappa_gen(T) = mu(T)*gamma/((gamma-1)Pr) = kappa(T)*gammaM2.
# Because mu/kappa now vary in space, div() below automatically includes the
# grad(mu), grad(kappa) contributions that horses3d's DG divergence also produces.
T_h    = gammaM2_sym * p / rho          # horses3d Temperature(Q)
tildeT = Tnorm_sym * T_h
suther = (1 + sdivTref_sym) / (tildeT + sdivTref_sym) * tildeT * sp.sqrt(tildeT)
mu_var    = eta_sym * suther
kappa_var = mu_var * gamma_sym / ((gamma_sym - 1) * Pr_sym)

# Viscous stress tensor components (variable viscosity)
tau_xx = mu_var*(2*ux - sp.Rational(2,3)*divU)
tau_yy = mu_var*(2*vy - sp.Rational(2,3)*divU)
tau_zz = mu_var*(2*wz - sp.Rational(2,3)*divU)
tau_xy = mu_var*(uy + vx)
tau_xz = mu_var*(uz + wx)
tau_yz = mu_var*(vz + wy)

# Temperature used in Fourier's law (ideal gas with R=1): T = p/rho
T = p / rho

# RHS expressions for continuity, momentum, energy
S_rho   = sp.diff(rho, t) + div(rho*u, rho*v, rho*w)
S_rhou  = sp.diff(rho*u, t) + div(rho*u*u + p, rho*u*v, rho*u*w) - div(tau_xx, tau_xy, tau_xz)
S_rhov  = sp.diff(rho*v, t) + div(rho*v*u, rho*v*v + p, rho*v*w) - div(tau_xy, tau_yy, tau_yz)
S_rhow  = sp.diff(rho*w, t) + div(rho*w*u, rho*w*v, rho*w*w + p) - div(tau_xz, tau_yz, tau_zz)
S_rhoe  = sp.diff(e_tot, t) \
          + div((e_tot + p)*u, (e_tot + p)*v, (e_tot + p)*w) \
          - div(tau_xx*u + tau_xy*v + tau_xz*w,
                tau_xy*u + tau_yy*v + tau_yz*w,
                tau_xz*u + tau_yz*v + tau_zz*w) \
          + div(-kappa_var*sp.diff(T,x), -kappa_var*sp.diff(T,y), -kappa_var*sp.diff(T,z))

rhs_exprs = [S_rho, S_rhou, S_rhov, S_rhow, S_rhoe]
var_names = ['IRHO', 'IRHOU', 'IRHOV', 'IRHOW', 'IRHOE']

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
    s = s.subs(gamma_sym, sp.Symbol('gamma'))
    s = s.subs(eta_sym, sp.Symbol('eta'))
    s = s.subs(Pr_sym, sp.Symbol('Pr'))
    s = s.subs(gammaM2_sym, sp.Symbol('gammaM2'))
    s = s.subs(sdivTref_sym, sp.Symbol('sdivTref'))
    s = s.subs(Tnorm_sym, sp.Symbol('Tnorm'))
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

    # Total energy (using gamma as a parameter)
    e_tot_expr = p_expr/(GAMMA-1) + sp.Rational(1,2)*rho_expr*(u_expr**2+v_expr**2+w_expr**2)

    # Convert to Fortran
    res_rho  = sp.fcode(rho_expr, assign_to='res(IRHO)',  source_format='free', standard=95)
    res_rhou = sp.fcode(rho_expr*u_expr, assign_to='res(IRHOU)', source_format='free', standard=95)
    res_rhov = sp.fcode(rho_expr*v_expr, assign_to='res(IRHOV)', source_format='free', standard=95)
    res_rhow = sp.fcode(rho_expr*w_expr, assign_to='res(IRHOW)', source_format='free', standard=95)
    res_rhoe = sp.fcode(e_tot_expr, assign_to='res(IRHOE)', source_format='free', standard=95)

    # Indent and combine
    def indent_block(code, ind='      '):
        return '\n'.join(ind + l for l in code.split('\n'))

    return f"""\
   subroutine evalManufacturedSolution(x, y, z, t, res)
!
!     Evaluates the manufactured solution at (x,y,z,t).
!     Conserved variables: rho, rhou, rhov, rhow, rhoe
!
!     RHO0 : reference density
!     GAMMA: specific heat ratio — must match control file.
!
      use SMConstants
      use PhysicsStorage_NS
      implicit none
      real(kind=RP), intent(in)  :: x, y, z, t
      real(kind=RP), intent(out) :: res(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP), parameter :: RHO0  = {RHO0}_RP
      real(kind=RP), parameter :: GAMMA = {GAMMA}_RP

{indent_block(res_rho)}
{indent_block(res_rhou)}
{indent_block(res_rhov)}
{indent_block(res_rhow)}
{indent_block(res_rhoe)}

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
!      ProblemFile.f90  —  HORSES3D Navier-Stokes (NS) MMS Convergence Study
!      Generated by generate_problemfile_NS.py
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

#ifdef NAVIERSTOKES
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
            real(kind=RP) :: t, rho_0, gamma, eta, Pr, gammaM2, sdivTref, Tnorm
            ! Sutherland's constant for air (horses3d default S_Sutherland = 198.6*5/9 K).
            real(kind=RP), parameter :: S_SUTHERLAND = 198.6_RP*5.0_RP/9.0_RP

#ifdef NAVIERSTOKES
            t     = time
            ! rho_0: reference density
            rho_0 = {RHO0}_RP
            gamma = thermodynamics_ % gamma
            eta   = dimensionless_  % mu        ! reference dynamic viscosity = 1/Re
            Pr    = dimensionless_  % Pr
            gammaM2 = dimensionless_ % gammaM2  ! gamma * Mach^2
            ! Sutherland's law (VariableConversion_NS.f90, get_laminar_mu_kappa):
            !   mu(T)    = eta * (1 + sdivTref)/(Tnorm*T_h + sdivTref) * (Tnorm*T_h)**1.5
            !   T_h      = gammaM2 * p / rho                (horses3d Temperature(Q))
            !   kappa(T) = mu(T) * gamma / ((gamma-1)*Pr)   [with Fourier T = p/rho]
            ! The reference temperature is taken at runtime from refValues_. The Sutherland
            ! constants live as module globals in PhysicsStorage_NS, but referencing them
            ! from this shared library would force linking non-PIC objects, so they are
            ! reconstructed here from horses3d's DEFAULTS (TRef_Sutherland = refValues%T,
            ! so Tnorm = 1). If you set "sutherland (reference) temperature" in the control
            ! file, update S_SUTHERLAND / the ratios below to match.
            Tnorm    = 1.0_RP                          ! refValues%T / TRef_Sutherland
            sdivTref = S_SUTHERLAND / refValues_ % T   ! S_Sutherland / TRef_Sutherland

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
            real(kind=RP)                      :: wi, wj, wk, wq
            real(kind=RP)                      :: error_mesh, vol
            !
            ! {NODE_TYPE} quadrature weights, orders 1..{MAX_ORDER}
            !
{weight_arrays_f90()}
            error_mesh = 0.0_RP
            vol        = 0.0_RP

#ifdef NAVIERSTOKES
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz ; DO j = 0, Ny ; DO i = 0, Nx
                  x   = mesh % elements(eID) % geom % x(:,i,j,k)
                  u_h = mesh % elements(eID) % storage % Q(:,i,j,k)
                  call evalManufacturedSolution(x(1), x(2), x(3), time, u_ex)

{weight_select_f90('Nx', 'i', 'wi')}
{weight_select_f90('Ny', 'j', 'wj')}
{weight_select_f90('Nz', 'k', 'wk')}
                  wq = wi*wj*wk * mesh % elements(eID) % geom % jacobian(i,j,k)
                  vol        = vol        + wq
                  error_mesh = error_mesh + wq * norm2(u_ex - u_h)**2
               END DO ; END DO ; END DO
            END DO

            ! Volume-averaged L2 (RMS) error: sqrt( integral(||err||^2) / volume ).
            ! The volume = integral(1) = sum of w*J is accumulated in the SAME loop,
            ! so this is correct for any mesh / domain (not just [-1,1]^3).
            error_mesh = sqrt(error_mesh / vol)

            ! --- Write L2 error to file for automated convergence study ---
            ! Format: nelems, P, NDOF, L2_error, t_final
            nelems = SIZE(mesh % elements)
            P      = mesh % elements(1) % Nxyz(1)
            NDOF   = nelems * (P+1)**3
            print *, "NAVIERSTOKES MMS L2 error: ", error_mesh
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

with open(OUTPUT_FILE, 'w') as f:
    f.write(pf)

print(f"Written {OUTPUT_FILE} ({len(pf.splitlines())} lines)", file=sys.stderr)
print(f"Next step: compile with 'make' in your HORSES3D setup directory.", file=sys.stderr)
