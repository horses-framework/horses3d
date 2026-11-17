#!/usr/bin/env python3
"""
Write <case>/SETUP/ProblemFile.f90 for an MMS case.

    python3 generate.py NS|Euler|NSSA|iNS|MU|CAA

Each case defines a manufactured solution, periodic on [-1,1]^3, and the
equations the solver discretises. The source term is the residual of the
manufactured solution, dq/dt + div F - S, emitted as it is, without symbolic
simplification. The physical parameters are read from the solver at run time,
except that NS and NSSA assume the default Sutherland constants, NSSA the
default Spalart-Allmaras constants, and MU per-fluid sound speeds
("fluid 1/2 sound speed square").
"""
import sys
from pathlib import Path
import sympy as sp

t, x, y, z = sp.symbols("t x y z")
pi = sp.pi
half = sp.Rational(1, 2)
SSS = sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)


def grad(f):
    return sp.diff(f, x), sp.diff(f, y), sp.diff(f, z)


def div(fx, fy, fz):
    return sp.diff(fx, x) + sp.diff(fy, y) + sp.diff(fz, z)


def lap(f):
    return div(*grad(f))


def cells(f, w=-2):
    """Cellular velocity field times f(t), divergence-free for w = -2."""
    return [sp.sin(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*f,
            sp.cos(pi*x)*sp.sin(pi*y)*sp.cos(pi*z)*f,
            w*sp.cos(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*f]


def stress(mu, V, bulk):
    """Viscous stress tensor mu (grad u + grad u^T) - bulk * mu * div(u) I."""
    g = [grad(a) for a in V]
    divu = g[0][0] + g[1][1] + g[2][2]
    return [[mu*(g[i][j] + g[j][i]) - (bulk*mu*divu if i == j else 0) for j in range(3)] for i in range(3)]


def momentum_flux(rho, V, p, tau):
    return [[rho*V[i]*a + (p if i == j else 0) - tau[i][j] for j, a in enumerate(V)] for i in range(3)]


def residual(q, F, S=None):
    """dq/dt + div F - S for each equation."""
    return [sp.diff(a, t) + div(*f) - s for a, f, s in zip(q, F, S or [0]*len(q))]


# Solver parameters, read at run time (Fortran expressions)
COMPRESSIBLE = {
    "gamma":   "thermodynamics_ % gamma",
    "eta":     "dimensionless_ % mu",
    "gammaM2": "dimensionless_ % gammaM2",
    "Pr":      "dimensionless_ % Pr",
    "Prt":     "dimensionless_ % Prt",
    # Default Sutherland constants (PhysicsStorage_NS): S = 110.4 K, TRef = refValues % T
    "s_suth":  "198.6_RP*5.0_RP/9.0_RP / refValues_ % T",
}


def compressible(rho, u, v, w, p, viscous=True, theta=None):
    """Euler, Navier-Stokes, or Navier-Stokes with Spalart-Allmaras (theta = nu_tilde)."""
    gamma, eta, gammaM2, Pr, Prt, s_suth = sp.symbols(list(COMPRESSIBLE))
    V = [u, v, w]
    E = p/(gamma - 1) + half*rho*(u**2 + v**2 + w**2)
    q = [rho, rho*u, rho*v, rho*w, E]
    names = ["IRHO", "IRHOU", "IRHOV", "IRHOW", "IRHOE"]
    tau, heat, S = [[0]*3]*3, [0]*3, None
    if viscous:
        # Sutherland's law at the solver's temperature T = gamma M^2 p/rho
        T = gammaM2*p/rho
        mu = (1 + s_suth)/(T + s_suth)*T**sp.Rational(3, 2)
        kappa = mu/((gamma - 1)*gammaM2/gamma*Pr)
        if theta is not None:
            # Spalart-Allmaras as in Turbulence.f90 (theta >= 0 branch). In a periodic box there
            # is no wall, the wall distance is sqrt(huge) and, of production and destruction,
            # only cb1 |omega| rho theta remains, so |omega| must not vanish anywhere.
            cv1, cb1, cb2, sigma = sp.Float(7.1), sp.Float(0.1355), sp.Float(0.622), sp.Rational(2, 3)
            g, gt = [grad(a) for a in V], grad(theta)
            omega = sp.sqrt((g[2][1] - g[1][2])**2 + (g[0][2] - g[2][0])**2 + (g[1][0] - g[0][1])**2)
            q.append(rho*theta)
            names.append("IRHOTHETA")
            theta_flux = [rho*theta*a - eta*(mu + rho*theta)/sigma*d for a, d in zip(V, gt)]
            S = [0]*5 + [cb1*omega*rho*theta + eta*cb2*rho*sum(d**2 for d in gt)/sigma]
            chi = rho*theta/mu
            mu_t = rho*theta*chi**3/(chi**3 + cv1**3)
            kappa += mu_t/((gamma - 1)*gammaM2/gamma*Prt)
            mu += mu_t
        tau = stress(eta*mu, V, sp.Rational(2, 3))
        heat = [-eta*kappa*d for d in grad(T)]
    F = [[rho*a for a in V], *momentum_flux(rho, V, p, tau),
         [(E + p)*V[j] - sum(tau[j][i]*V[i] for i in range(3)) + heat[j] for j in range(3)]]
    if theta is not None:
        F.append(theta_flux)
    return q, residual(q, F, S), names, COMPRESSIBLE


# Smooth field at full amplitude for Euler and NSSA. The velocity has the vorticity of
# (sin(pi z), cos(pi z), 0) plus a smaller 3D part, so |omega| > 0.
CELLS = cells(1)
SUSTAINED = dict(rho=1 + 0.2*sp.sin(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.cos(t),
                 u=0.3*sp.cos(t)*(sp.sin(pi*z) + 0.1*CELLS[0]),
                 v=0.3*sp.cos(t)*(sp.cos(pi*z) + 0.1*CELLS[1]),
                 w=0.3*sp.cos(t)*(0.1*CELLS[2]),
                 p=1 + 0.1*SSS*sp.cos(t))


def navier_stokes():
    return compressible(sp.Integer(1), *cells(sp.sin(t)), 1 + 0.1*SSS*sp.cos(t))


def euler():
    return compressible(**SUSTAINED, viscous=False)


def spalart_allmaras():
    theta = 1 + 0.5*sp.cos(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.cos(t)
    return compressible(**SUSTAINED, theta=theta)


def incompressible():
    """Artificial compressibility: dp/dt + rho0 c0^2 div(u) = 0, with rho = 1."""
    mu, c0_sq = sp.symbols("mu c0_sq")
    V, p = cells(sp.sin(t)), 1 + 0.1*SSS*sp.cos(t)
    q = [sp.Integer(1), *V, p]
    F = [V, *momentum_flux(1, V, p, stress(mu, V, 0)), [c0_sq*a for a in V]]
    params = {"mu": "dimensionless_ % mu(1)", "c0_sq": "thermodynamics_ % rho0c02"}
    return q, residual(q, F), ["INSRHO", "INSRHOU", "INSRHOV", "INSRHOW", "INSP"], params


def multiphase():
    """Cahn-Hilliard + incompressible NS with artificial compressibility, two fluids."""
    rho_1, rho_2, eta_1, eta_2, cs_1, cs_2, M0, sigma, eps = sp.symbols("rho_1 rho_2 eta_1 eta_2 cs_1 cs_2 M0 sigma eps")
    # cos(t) so the field starts at full amplitude: the explicit time step is
    # limited by Cahn-Hilliard (~h^4) and the runs only reach t ~ 1e-4.
    c = half*(1 + sp.cos(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.cos(t))
    V, p = cells(sp.cos(t)), SSS*sp.sin(t)
    # Mixture density, viscosity and sound speed as in the solver
    rho = c*rho_1 + (1 - c)*rho_2
    sqrho = sp.sqrt(rho)
    eta = c*eta_1 + (1 - c)*eta_2
    cs = c*cs_1 + (1 - c)*cs_2
    # Chemical potential
    mu_c = (sigma/eps)*(12*c**2*(2*c - 2) + 24*c*(1 - c)**2) - sp.Rational(3, 2)*sigma*eps*lap(c)
    tau = stress(eta, V, 0)
    # Not in conservative form: the solver splits the convection half in divergence and
    # half in advective form (Physics_MU.f90), and the variables are c, sqrt(rho) u, p.
    source = [sp.diff(c, t) + div(*[c*a for a in V]) - M0*lap(mu_c)]
    source += [sqrho*sp.diff(sqrho*V[i], t) + half*div(*[rho*a*V[i] for a in V])
               + half*rho*sum(a*d for a, d in zip(V, grad(V[i]))) + c*grad(mu_c)[i] + grad(p)[i] - div(*tau[i])
               for i in range(3)]
    source += [sp.diff(p, t) + cs**2*rho*div(*V)]
    params = {"rho_1": "dimensionless_ % rho(1)", "rho_2": "dimensionless_ % rho(2)",
              "eta_1": "dimensionless_ % mu(1)", "eta_2": "dimensionless_ % mu(2)",
              "cs_1": "sqrt(dimensionless_ % invMa2(1) / dimensionless_ % rho(1))",
              "cs_2": "sqrt(dimensionless_ % invMa2(2) / dimensionless_ % rho(2))",
              "M0": "multiphase_ % M0", "sigma": "multiphase_ % sigma", "eps": "multiphase_ % eps"}
    return [c, *[sqrho*a for a in V], p], source, ["IMC", "IMSQRHOU", "IMSQRHOV", "IMSQRHOW", "IMP"], params


def acoustics():
    """APE-4 perturbations on the uniform base flow of the control file (no Lamb vector)."""
    rho0, u0, v0, w0, c0_sq = sp.symbols("rho0 u0 v0 w0 c0_sq")
    V, p = cells(sp.cos(t), w=1), SSS*sp.sin(t)
    phi = u0*V[0] + v0*V[1] + w0*V[2] + p/rho0
    q = [sp.Integer(0), *V, p]
    F = [[0, 0, 0], [phi, 0, 0], [0, phi, 0], [0, 0, phi], [rho0*c0_sq*a + b*p for a, b in zip(V, [u0, v0, w0])]]
    params = {"rho0": "Qbase(IBRHO)", "u0": "Qbase(IBU)", "v0": "Qbase(IBV)", "w0": "Qbase(IBW)", "c0_sq": "Qbase(IBA2)"}
    return q, residual(q, F), ["ICAARHO", "ICAAU", "ICAAV", "ICAAW", "ICAAP"], params


# Equations, solver, and the names of the variables in the error file
CASES = {
    "NS":    (navier_stokes,    "ns",   ["rho", "rhou", "rhov", "rhow", "rhoe"]),
    "Euler": (euler,            "ns",   ["rho", "rhou", "rhov", "rhow", "rhoe"]),
    "NSSA":  (spalart_allmaras, "nssa", ["rho", "rhou", "rhov", "rhow", "rhoe", "rhotheta"]),
    "iNS":   (incompressible,   "ins",  ["rho", "rhou", "rhov", "rhow", "p"]),
    "MU":    (multiphase,       "mu",   ["c", "sqrhou", "sqrhov", "sqrhow", "p"]),
    "CAA":   (acoustics,        "caa",  ["rho", "u", "v", "w", "p"]),
}
MACROS = {"ns": "NAVIERSTOKES", "nssa": "SPALARTALMARAS", "ins": "INCNS", "mu": "MULTIPHASE", "caa": "ACOUSTIC"}


def routine(array, names, exprs, params):
    """Body of MMS_State or MMS_Source: the parameters used, then one assignment per variable."""
    exprs = [sp.sympify(e) for e in exprs]
    used = sorted({s.name for e in exprs for s in e.free_symbols} & params.keys())
    lines = [f"real(kind=RP) :: {', '.join(used)}"] if used else []
    lines += [f"{name} = {params[name]}" for name in used]
    fortran_names = {pi: sp.Symbol("PI"), x: sp.Symbol("x(1)"), y: sp.Symbol("x(2)"), z: sp.Symbol("x(3)")}
    for name, e in zip(names, exprs):
        lines += sp.fcode(e.subs(fortran_names), assign_to=f"{array}({name})", source_format="free",
                          standard=95, contract=False).splitlines()
    return "\n".join("         " + l for l in lines)


def generate(case):
    equations, solver, _ = CASES[case]
    exact, source, names, params = equations()
    here = Path(__file__).resolve().parent
    pf = (here / "ProblemFile.f90.in").read_text()
    pf = pf.replace("@SOLVER@", case).replace("@MACRO@", MACROS[solver])
    pf = pf.replace("@STATE@", routine("Q", names, exact, params))
    pf = pf.replace("@SOURCE@", routine("S", names, source, params))
    out = here / case / "SETUP" / "ProblemFile.f90"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(pf)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in CASES:
        sys.exit(f"usage: {sys.argv[0]} {'|'.join(CASES)}")
    print(f"Written {generate(sys.argv[1])}")
