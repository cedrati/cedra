"""
Cedra constant: core definitions and analysis helpers.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Only the Python standard library is used. Constants are computed with
`decimal` at 60 significant digits so that {n * alpha} stays exact to
far beyond double precision (plain floats drift by ~1e-7 at n = 1e9).
"""

from decimal import Decimal, ROUND_FLOOR, getcontext
import math

getcontext().prec = 60

_SQRT2 = Decimal(2).sqrt()
_SQRT3 = Decimal(3).sqrt()
_SQRT5 = Decimal(5).sqrt()

# C = sqrt(3) + sqrt(2) + sqrt(1/2) - 2 = sqrt(3) + (3/2) sqrt(2) - 2
CEDRA_DEC = _SQRT3 + _SQRT2 + (Decimal(1) / 2).sqrt() - 2
ALPHA_DEC = CEDRA_DEC.ln()          # rotation number alpha = ln C (transcendental)
INV_PHI_DEC = (_SQRT5 - 1) / 2      # 1/phi
SQRT2_MINUS_1_DEC = _SQRT2 - 1

CEDRA = float(CEDRA_DEC)
ALPHA = float(ALPHA_DEC)
INV_PHI = float(INV_PHI_DEC)
PHI = 1 + INV_PHI

# Minimal polynomial of C over Q (integer coefficients, highest degree first)
MINIMAL_POLYNOMIAL = (4, 32, 36, -112, -167)


def eval_poly(coeffs, x):
    """Horner evaluation, works with Decimal or float."""
    acc = 0
    for c in coeffs:
        acc = acc * x + c
    return acc


def frac(x):
    """Fractional part {x} = x - floor(x), valid for negative x too."""
    if isinstance(x, Decimal):
        return x - x.to_integral_value(rounding=ROUND_FLOOR)
    return x - math.floor(x)


def kronecker(n_points, alpha=ALPHA_DEC, start=1):
    """Kronecker (irrational rotation) sequence {n * alpha}, n = start..start+n_points-1."""
    alpha = Decimal(alpha)
    return [float(frac(n * alpha)) for n in range(start, start + n_points)]


def window_word(n_points, window=INV_PHI_DEC, alpha=ALPHA_DEC, start=1):
    """Binary coding of the rotation: 1 if {n * alpha} < window, else 0."""
    alpha, window = Decimal(alpha), Decimal(window)
    return [1 if frac(n * alpha) < window else 0
            for n in range(start, start + n_points)]


def model_set(n_points, window=INV_PHI_DEC, alpha=ALPHA_DEC, start=1):
    """Cut-and-project (model) set S = {n | {n * alpha} < window} restricted to a range."""
    word = window_word(n_points, window, alpha, start)
    return [start + i for i, bit in enumerate(word) if bit]


def continued_fraction(x, n_terms=20):
    """First partial quotients of a positive Decimal."""
    x = Decimal(x)
    terms = []
    for _ in range(n_terms):
        q = int(x.to_integral_value(rounding=ROUND_FLOOR))
        terms.append(q)
        x -= q
        if x == 0:
            break
        x = 1 / x
    return terms


def convergent_denominators(x, n_terms=12):
    """Denominators q_k of the continued fraction convergents p_k / q_k of x."""
    q_prev, q = 0, 1
    dens = []
    for a in continued_fraction(x, n_terms):
        if dens:
            q_prev, q = q, a * q + q_prev
        dens.append(q)
    return dens


def distance_to_integer(x):
    """||x||: distance from x to the nearest integer."""
    f = frac(x)
    return min(f, 1 - f)


def factor_complexity(word, k):
    """Number of distinct factors (sub-words) of length k."""
    s = ''.join(map(str, word))
    return len({s[i:i + k] for i in range(len(s) - k + 1)})


def smallest_period(word, max_period=None):
    """Smallest p such that word[i] == word[i + p] for all i, or None."""
    n = len(word)
    max_period = max_period or n // 2
    for p in range(1, max_period + 1):
        if all(word[i] == word[i + p] for i in range(n - p)):
            return p
    return None


def star_discrepancy(points):
    """Exact 1D star discrepancy D*_N of points in [0, 1)."""
    s = sorted(points)
    n = len(s)
    return max(max((i + 1) / n - x, x - i / n) for i, x in enumerate(s))


def chi_squared(points, bins=20):
    """Pearson chi-squared statistic against the uniform law on [0, 1)."""
    counts = [0] * bins
    for x in points:
        counts[min(int(x * bins), bins - 1)] += 1
    expected = len(points) / bins
    return sum((c - expected) ** 2 / expected for c in counts)


def serial_correlation(data, lag=1):
    """Pearson correlation between data[i] and data[i + lag]."""
    x, y = data[:-lag], data[lag:]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    vx = sum((a - mx) ** 2 for a in x)
    vy = sum((b - my) ** 2 for b in y)
    return cov / math.sqrt(vx * vy)


def rotation_autocorrelation(lag, alpha=ALPHA):
    """Exact asymptotic autocorrelation of {n * alpha} at a given lag: 1 - 6 f (1 - f), f = {lag * alpha}."""
    f = frac(lag * float(alpha))
    return 1 - 6 * f * (1 - f)
