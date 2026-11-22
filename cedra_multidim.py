"""
Cedra in several dimensions: numerical check of docs/higher-dimensions.md.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

The four Galois conjugates C_i = +-sqrt(3) +- (3/2) sqrt(2) - 2 give
1, ln|C_1|, ..., ln|C_4| linearly independent over Q (Theorem D), so
(n ln|C_1|, ..., n ln|C_4|) is equidistributed in the 4-torus. This
script checks the number theory behind the proof and measures how good
the resulting point sets are compared with standard constructions.
"""

from decimal import Decimal
import random

from cedra import (CONJUGATES_DEC, LOG_CONJUGATES_DEC, MINIMAL_POLYNOMIAL, eval_poly,
                   halton, kronecker_nd, l2_star_discrepancy)

PLASTIC = Decimal("1.324717957244746025960908854478097340734404056901733364534")
SEED = 1


def number_theory():
    ok = True
    print("Galois conjugates of C and their logarithms")
    for i, (c, l) in enumerate(zip(CONJUGATES_DEC, LOG_CONJUGATES_DEC), start=1):
        print(f"  C{i} = {float(c):+.12f}   ln|C{i}| = {float(l):.12f}")

    norm = 16 * CONJUGATES_DEC[0] * CONJUGATES_DEC[1] * CONJUGATES_DEC[2] * CONJUGATES_DEC[3]
    good = abs(norm + 668) < Decimal("1e-40")
    ok &= good
    print(f"\n[{'PASS' if good else 'FAIL'}] N(2C) = 16 * C1 C2 C3 C4 = {round(norm)} = -4 * 167")

    p = 167
    squares = pow(2, (p - 1) // 2, p) == 1 and pow(3, (p - 1) // 2, p) == 1
    roots = [x for x in range(p) if eval_poly(MINIMAL_POLYNOMIAL, x) % p == 0]
    good = squares and len(roots) == 4
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] 2 and 3 are squares mod 167, P has 4 roots mod 167: {roots}")
    print("  -> 167 splits completely in Q(sqrt 2, sqrt 3); each C_i is divisible by its own prime above 167")
    return ok


def discrepancy_tables():
    rng = random.Random(SEED)
    lc = LOG_CONJUGATES_DEC

    print("\nL2 star discrepancy in 2D (lower is better)")
    sets_2d = {
        "Cedra (ln C1, ln|C2|)": lambda n: kronecker_nd(n, (lc[0], lc[1])),
        "Cedra (ln C1, ln|C3|)": lambda n: kronecker_nd(n, (lc[0], lc[2])),
        "R2 (plastic number)": lambda n: kronecker_nd(n, (1 / PLASTIC, 1 / PLASTIC ** 2)),
        "(sqrt 2, sqrt 3)": lambda n: kronecker_nd(n, (Decimal(2).sqrt(), Decimal(3).sqrt())),
        "Halton (2, 3)": lambda n: halton(n, (2, 3)),
        "pseudo-random": lambda n: [(rng.random(), rng.random()) for _ in range(n)],
    }
    sizes = (500, 2000)
    print(f"  {'point set':<24}" + "".join(f"{'N=' + str(n):>12}" for n in sizes))
    for name, gen in sets_2d.items():
        print(f"  {name:<24}" + "".join(f"{l2_star_discrepancy(gen(n)):>12.3e}" for n in sizes))

    n = 1500
    print(f"\nL2 star discrepancy in 4D, N = {n}")
    sets_4d = {
        "Cedra (all 4 conjugates)": kronecker_nd(n, lc),
        "(sqrt 2, 3, 5, 7)": kronecker_nd(n, [Decimal(k).sqrt() for k in (2, 3, 5, 7)]),
        "Halton (2, 3, 5, 7)": halton(n, (2, 3, 5, 7)),
        "pseudo-random": [tuple(rng.random() for _ in range(4)) for _ in range(n)],
    }
    for name, pts in sets_4d.items():
        print(f"  {name:<24}{l2_star_discrepancy(pts):>12.3e}")


def simultaneous_records(a, b, n_max):
    """Records of max(||n a||, ||n b||) and the normalised value n * max^2."""
    a, b = float(a), float(b)
    best, records = 1.0, []
    for n in range(1, n_max + 1):
        x, y = n * a, n * b
        m = max(abs(x - round(x)), abs(y - round(y)))
        if m < best:
            best = m
            records.append((n, n * m * m))
    return records


def approximation_tables(n_max=2_000_000):
    print(f"\nSimultaneous approximation up to n = {n_max:,}: records of n * max(||n a||, ||n b||)^2")
    print("  (bounded below for a badly approximable pair; Dirichlet guarantees values <= 1)")
    lc = LOG_CONJUGATES_DEC
    pairs = {
        "Cedra (ln C1, ln|C2|)": (lc[0], lc[1]),
        "R2 (plastic number)": (1 / PLASTIC, 1 / PLASTIC ** 2),
        "(sqrt 2, sqrt 3)": (Decimal(2).sqrt(), Decimal(3).sqrt()),
    }
    for name, (a, b) in pairs.items():
        rec = [r for r in simultaneous_records(a, b, n_max) if r[0] > 100]
        worst = min(rec, key=lambda r: r[1])
        print(f"  {name:<24} smallest {worst[1]:.4f} at n = {worst[0]:>8,}   "
              f"last records: {', '.join(f'{n:,}' for n, _ in rec[-4:])}")


if __name__ == "__main__":
    ok = number_theory()
    discrepancy_tables()
    approximation_tables()
    print("\nNumber theory checks passed." if ok else "\nSome checks FAILED.")
