#!/usr/bin/env python3
"""Exact rational checks for the paper "A lower bound for ropelength in terms of
crossing number". Python standard library only.

This is the computational part of an analytic proof, NOT a formalization of
geometric measure theory. The geometric lemmas are proved in the paper.
No float, numerical quadrature, sampled maximum, or downloaded data determines
any PASS. perf_counter is used only to report elapsed time.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json
from pathlib import Path
from time import perf_counter


def require(condition, label):
    if not condition:
        raise RuntimeError("FAILED: " + label)


def atan_interval(x, n):
    """Alternating series with n terms and rigorous next-term enclosure."""
    require(0 < x < 1 and n > 0, "atan series domain")
    s = sum(((-1) ** j * x ** (2*j+1) / (2*j+1)
             for j in range(n)), Q(0))
    next_term = (-1) ** n * x ** (2*n+1) / (2*n+1)
    return min(s, s + next_term), max(s, s + next_term)


def log4_interval(n):
    """log(4)=4 atanh(1/3); bound the positive series tail geometrically."""
    z = Q(1, 3)
    s = 4 * sum((z ** (2*j+1) / (2*j+1) for j in range(n)), Q(0))
    tail = 4*z ** (2*n+1) / ((2*n+1)*(1-z*z))
    return s, s + tail


def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def ladd(*polys):
    """Add Laurent polynomials (integer exponent -> rational coefficient)."""
    result = {}
    for p in polys:
        for e, c in p.items():
            result[e] = result.get(e, Q(0)) + c
    return {e: c for e, c in result.items() if c}


def deriv(p):
    return {e-1: e*c for e, c in p.items() if e and c}


def peval(p, x):
    return sum((c*x**e for e, c in p.items()), Q(0))


def approx(q):
    """Presentation only; all decisions above/below use Fraction."""
    with localcontext() as ctx:
        ctx.prec = 32
        return str(Decimal(q.numerator)/Decimal(q.denominator))


def main():
    start = perf_counter()
    k = Q(50657, 10000)
    u = Q(138629436112, 100000000000)
    claim = Q(531, 250)  # 2.124 exactly.

    # Machin: 4 atan(1/5)-atan(1/239)=pi/4. Exact tangent check;
    # the angle lies in (0,1), so there is no tangent branch ambiguity.
    z = (Q(1), Q(0))
    for _ in range(4):
        z = cmul(z, (Q(1), Q(1, 5)))
    z = cmul(z, (Q(1), Q(-1, 239)))
    require(z[0] == z[1] and z[0] > 0, "Machin tangent identity")
    a_lo, a_hi = atan_interval(Q(1, 5), 16)
    b_lo, b_hi = atan_interval(Q(1, 239), 8)
    require(0 < 4*a_lo-b_hi <= 4*a_hi-b_lo < 1,
            "Machin angle branch (0,1)")
    pi_lo, pi_hi = 16*a_lo-4*b_hi, 16*a_hi-4*b_lo
    require(3 < pi_lo < pi_hi < Q(22, 7), "pi enclosure")

    log_lo, log_hi = log4_interval(16)
    require(log_lo < log_hi < u, "log(4) certified upper bound")

    # Integrate (8/3)*(t+1)^3/t^3. Its primitive is
    # J(t)=(8/3)t+8 log(t)-8/t-4/(3t^2).
    j = {1: Q(8, 3), -1: Q(-8), -2: Q(-4, 3)}
    integrand = {0: Q(8, 3), -1: Q(8), -2: Q(8), -3: Q(8, 3)}
    require(ladd(deriv(j), {-1: Q(8)}) == integrand,
            "radial integral primitive")
    # The tail ell/R^2 for ell=(4/3)(R+1)^3.
    tail = {1: Q(4, 3), 0: Q(4), -1: Q(4), -2: Q(4, 3)}
    f = ladd(j, {0: -peval(j, Q(2))}, tail)
    require(f == {1: Q(4), 0: Q(3), -1: Q(-4)},
            "F(R)=4R+8 log(R/2)+3-4/R")
    require(peval(f, Q(2)) == 9, "continuity at ell=36, R=2")

    # Concavity: log(R/2) <= log(4)+(R-8)/8.
    # This replaces F by 5R+8u-5-4/R.
    majorant = ladd(f, {0: 8*u-8, 1: Q(1)})
    require(majorant == {1: Q(5), 0: 8*u-5, -1: Q(-4)},
            "logarithmic tangent majorant")
    difference = ladd({1: k, 0: k}, {e: -c for e, c in majorant.items()})
    p = {e+1: c for e, c in difference.items()}
    aa, bb, cc = k-5, k+5-8*u, Q(4)
    require(p == {2: aa, 1: bb, 0: cc}, "quadratic reduction")
    gap = 4*aa*cc-bb*bb
    require(aa > 0 and gap > 0, "global positive quadratic")
    # 4a P(R) = (2aR+b)^2 + gap, for EVERY real R.
    require(4*aa*cc == bb*bb+gap, "completed-square certificate")

    # Small-length case: (ell/4)^3 <= (k*(3/4)^(1/3)*ell^(1/3))^3.
    require(Q(36)**2 < 48*k**3, "all 0<ell<=36")

    # If A=k*(3/4)^(1/3), C<=A ell^(4/3)/(4pi).
    # The corresponding lower bound B obeys B^4=256*pi^3/(3*k^3).
    bound4_lo = 256*pi_lo**3/(3*k**3)
    bound4_hi = 256*pi_hi**3/(3*k**3)
    require(claim**4 < bound4_lo, "C_22a >= 531/250")
    # Enclose the unrounded exact-expression bound as well.
    b_lo = Q(212403897858, 100000000000)
    b_hi = Q(212403897859, 100000000000)
    require(b_lo**4 < bound4_lo and bound4_hi < b_hi**4,
            "unrounded bound interval")
    for denominator in (Q(11), Q(1067, 100)):
        require(claim**4 > (4*pi_hi/denominator)**3,
                "strict improvement over denominator " + str(denominator))
    require(claim**4 > (4*pi_hi/Q(646, 100))**3,
            "strict improvement over Diao-Ernst Theorem 3.2 (knots)")

    # Negative control: the tempting rounded coefficient 4.6 is invalid.
    # At R=8 its cube must exceed F(8)^3/ell, but does not.
    f8_lower = peval(f, Q(8)) + 8*log_lo
    ell8 = Q(4, 3)*9**3
    require(f8_lower**3 > Q(23, 5)**3*ell8,
            "reject unsafe coefficient A=4.6")

    certificate = {
        "verified": True,
        "statement": "C_22a >= 531/250 = 2.124 for knots and links (unit-radius convention).",
        "exact_bound": "(256*pi^3/(3*(50657/10000)^3))^(1/4)",
        "exact_bound_decimal_enclosure": [approx(b_lo), approx(b_hi)],
        "arithmetic": "Python fractions.Fraction (exact integer/rational arithmetic)",
        "coverage": "All ell>0 and all R>=2 via analytic reduction and global quadratic certificate; no sampling.",
        "scope": "Rational inequalities of the paper; the geometric lemmas and the integration argument are proved in the paper, not checked here.",
        "k": str(k),
        "log4_upper": str(u),
        "log4_upper_margin": str(u-log_hi),
        "quadratic_gap": str(gap),
        "pi_interval": [str(pi_lo), str(pi_hi)],
        "claimed_bound_fourth_power_margin": str(bound4_lo-claim**4),
        "earlier_denominators_checked": ["11", "10.67"],
        "additional_literature_comparison": "Diao-Ernst Theorem 3.2: denominator 6.46 (knots)",
        "negative_control": "A=4.6 rejected at R=8",
        "elapsed_seconds": round(perf_counter()-start, 6),
    }
    out = Path(__file__).resolve().parent / "certificate.json"
    out.write_text(json.dumps(certificate, ensure_ascii=False, indent=2) + "\n")
    print("PASS: C_22a >= 531/250 = 2.124")
    print("Exact-expression lower bound in [" + approx(b_lo) + ", " + approx(b_hi) + "]")
    print("All rational checks passed; the geometric part is proved in the paper.")
    print("Verification seconds:", certificate["elapsed_seconds"])


if __name__ == "__main__":
    main()
