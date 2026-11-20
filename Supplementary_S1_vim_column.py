"""
Supplementary Material S1 / Қосымша материал S1 / Дополнительный материал S1

KK: Сыни күш: өз салмағы бар консольді сырық және Винклер іргетасындағы сырық.
RU: Критическая сила: консольный стержень с собственным весом и стержень на основании Винклера.
EN: Critical load: self-weighted cantilever column and column on a Winkler foundation.

Reference high-accuracy numerical solution: shooting method, scipy.integrate.solve_ivp,
tolerances rtol, atol within 1e-10 .. 1e-12 (DOP853 / RK45).
VIM roots: mpmath, 40 digits. Symbolic derivations: SymPy.
Run: python vim_column_S1.py   (about 1-2 minutes)
"""
import sys
import numpy as np
import scipy
import sympy as sp
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

mp.mp.dps = 40

# ---------------------------------------------------------------
# 1. SymPy: VIM for the self-weight problem
#    theta'' + (lam + mu(1-xi)) theta = 0, theta(0)=0, theta'(1)=0
# ---------------------------------------------------------------
def derive_selfweight_vim():
    xi, s, lam, mu = sp.symbols('xi s lam mu', positive=True)
    w = sp.sqrt(lam)
    Lam = sp.sin(w * (s - xi)) / w            # Lagrange multiplier

    def step(th):
        ths = th.subs(xi, s)
        integrand = Lam * (sp.diff(ths, s, 2) + (lam + mu * (1 - s)) * ths)
        return sp.simplify(th + sp.integrate(sp.expand(integrand), (s, 0, xi)))

    th1 = step(xi)                             # theta_0 = A*xi, A=1
    th2 = step(th1)
    f1 = sp.simplify(sp.diff(th1, xi).subs(xi, 1))
    f2 = sp.simplify(sp.diff(th2, xi).subs(xi, 1))
    return th1, f1, f2, (xi, lam, mu)

# Closed forms of the VIM eigenvalue equations (for fast root finding)
def F1(l, m):                                  # VIM-1: theta_1'(1)=0
    l = mp.mpc(l); w = mp.sqrt(l)
    return mp.re(mp.cos(w) + m*(1 + mp.cos(w))/l - 2*m*mp.sin(w)/l**1.5)

def F2(l, m):                                  # VIM-2: theta_2'(1)=0
    l = mp.mpc(l); w = mp.sqrt(l); c = mp.cos(w); sn = mp.sin(w)
    return mp.re(c - m/4*(c/l + sn/w - sn/(l*w))
                 + m**2*(-7*c/(4*l**2) + 8*(c-1)/l**3
                         + 23*sn/(4*l**2*w) - sn/(4*l*w)))

# ---------------------------------------------------------------
# 2. Reference numerical solution (shooting method)
# ---------------------------------------------------------------
def shoot_selfweight(lam, mu):
    sol = solve_ivp(lambda x, u: [u[1], -(lam + mu*(1-x))*u[0]],
                    [0, 1], [0, 1], rtol=1e-12, atol=1e-12, method='DOP853')
    return sol.y[1, -1]

def shoot_winkler_det(lam, mu, kap):
    def rhs(x, u):
        Y, Y1, Y2, Y3 = u
        return [Y1, Y2, Y3, -(lam + mu*(1-x))*Y2 + mu*Y1 - kap*Y]
    rows = []
    for ic in ([0, 0, 1, 0], [0, 0, 0, 1]):
        s = solve_ivp(rhs, [0, 1], ic, rtol=1e-10, atol=1e-12)
        Y, Y1, Y2, Y3 = s.y[:, -1]
        rows.append([Y2, Y3 + lam*Y1])
    return np.linalg.det(np.array(rows))

def first_root_np(fun, lo, hi, n=300):
    xs = np.linspace(lo, hi, n); v = [fun(x) for x in xs]
    for i in range(n - 1):
        if v[i]*v[i+1] < 0:
            return brentq(fun, xs[i], xs[i+1], xtol=1e-13)

def first_root_mp(F, m, lo, hi, n=300):
    xs = [mp.mpf(lo) + (mp.mpf(hi) - lo)*i/n for i in range(n + 1)]
    v = [F(x, m) for x in xs]
    for i in range(n):
        if v[i]*v[i+1] < 0:
            return float(mp.findroot(lambda x: F(x, m), (xs[i], xs[i+1]),
                                     solver='anderson', tol=1e-25))

# ---------------------------------------------------------------
# 3. SymPy: VIM-1 for the Winkler problem -> determinant a11*a22 - a12*a21
#    Y'''' + (lam + mu(1-xi))Y'' - mu Y' + kap Y = 0
#    Y(0)=Y'(0)=0, Y''(1)=0, Y'''(1)+lam Y'(1)=0
# ---------------------------------------------------------------
def derive_winkler_vim1():
    xi, s, lam, mu, kap, A, B = sp.symbols('xi s lam mu kap A B', positive=True)
    w = sp.sqrt(lam)
    Lam = (s - xi)/lam - sp.sin(w*(s - xi))/(lam*w)
    Y0 = A*xi**2 + B*xi**3
    Y0s = Y0.subs(xi, s)
    res = (sp.diff(Y0s, s, 4) + (lam + mu*(1 - s))*sp.diff(Y0s, s, 2)
           - mu*sp.diff(Y0s, s) + kap*Y0s)
    Y1 = Y0 + sp.integrate(sp.expand(Lam*res), (s, 0, xi))
    c1 = sp.diff(Y1, xi, 2).subs(xi, 1)
    c2 = (sp.diff(Y1, xi, 3) + lam*sp.diff(Y1, xi)).subs(xi, 1)
    a11, a12 = sp.simplify(sp.diff(c1, A)), sp.simplify(sp.diff(c1, B))
    a21, a22 = sp.simplify(sp.diff(c2, A)), sp.simplify(sp.diff(c2, B))
    return (a11, a12, a21, a22), (lam, mu, kap)

# ---------------------------------------------------------------
if __name__ == "__main__":
    print("python", sys.version.split()[0], "| numpy", np.__version__,
          "| scipy", scipy.__version__, "| sympy", sp.__version__, "| mpmath", mp.__version__)
    mu_star = brentq(lambda m: shoot_selfweight(0.0, m), 5, 10, xtol=1e-13)
    print(f"Greenhill mu* (reference) = {mu_star:.4f}")
    lam_s, mu_s = sp.symbols('lam mu')
    f1s = sp.cos(sp.sqrt(lam_s)) + mu_s*(1+sp.cos(sp.sqrt(lam_s)))/lam_s \
          - 2*mu_s*sp.sin(sp.sqrt(lam_s))/lam_s**sp.Rational(3, 2)
    print("lim lam->0, VIM-1:", sp.limit(f1s, lam_s, 0), "-> mu*_1 = 6")
    print("VIM-2: 1 - mu/6 + mu^2/180 = 0 -> mu*_2 =",
          round(float(sp.solve(1 - mu_s/6 + mu_s**2/180, mu_s)[0]), 4))

    print("\nTable 1 (mu < mu*)")
    print(f"{'mu':>4} {'ref':>9} {'VIM-1':>9} {'err1%':>8} {'VIM-2':>9} {'err2%':>8}")
    for m in [0, 1, 2, 3, 4, 5, 6, 7]:
        ex = first_root_np(lambda l: shoot_selfweight(l, m), 0.05, 30)
        v1 = first_root_mp(F1, m, 0.002, 10)   # window <10 excludes the 2nd eigenvalue
        v2 = first_root_mp(F2, m, 0.002, 10)
        f = lambda v: f"{v:9.4f}" if v is not None else f"{'(<=0)':>9}"
        e = lambda v: f"{(v/ex-1)*100:8.2f}" if v is not None else f"{'-':>8}"
        print(f"{m:4d} {f(ex)} {f(v1)} {e(v1)} {f(v2)} {e(v2)}")

    (a11, a12, a21, a22), (lam, mu, kap) = derive_winkler_vim1()
    print("\na21 =", a21, "| a22 =", a22)
    D = sp.lambdify((lam, mu, kap), a11*a22 - a12*a21, 'mpmath')
    print("\nTable 3 (Winkler)")
    print(f"{'mu':>4} {'kappa':>6} {'ref':>9} {'VIM-1':>9} {'err%':>7}")
    for m, k in [(0, 0), (0, 10), (0, 27.78), (0, 50), (1, 0), (1, 10), (1, 27.78), (2, 20)]:
        ex = first_root_np(lambda l: shoot_winkler_det(l, m, k), 1e-3, 20)
        v = first_root_mp(lambda l, _m: D(l, mp.mpf(m), mp.mpf(k)), 0, 1.0, 25, 400)
        print(f"{m:4} {k:6} {ex:9.4f} {v:9.4f} {(v/ex-1)*100:7.2f}")

    E, b, h, l = 2.1e11, 0.1, 0.1, 3.0
    A_, I = b*h, b*h**3/12
    EI = E*I; q = 7850*A_*9.81; P = 350e3; kw = 6e6*b
    lam_p, mu_p, kap_p = P*l**2/EI, q*l**3/EI, kw*l**4/EI
    print(f"\nExample: lam={lam_p:.3f}, mu={mu_p:.5f}, kappa={kap_p:.2f}")
    ex = first_root_np(lambda x: shoot_winkler_det(x, mu_p, kap_p), 1e-3, 20)
    v1 = first_root_mp(lambda x, _m: D(x, mp.mpf(mu_p), mp.mpf(kap_p)), 0, 1.0, 25, 400)
    k_ = EI/l**2/1e3
    print(f"lam_cr ref={ex:.3f} -> Pcr={ex*k_:.1f} kN (n={ex*k_*1e3/P:.2f});"
          f" VIM-1={v1:.3f} -> {v1*k_:.1f} kN")
    print(f"Euler (cantilever): {np.pi**2/4*k_:.1f} kN")
