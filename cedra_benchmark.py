"""
Cedra benchmark: reproducible comparison of Kronecker sequences.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Compares {n alpha} for alpha = ln C, 1/phi and sqrt(2) - 1, plus a
pseudo-random generator. Lower star discrepancy = better for quasi-Monte
Carlo. A very low chi-squared means "more regular than random", not
"more random": a true uniform generator gives chi^2 ~ 19 with 20 bins.
"""

import random

from cedra import (ALPHA, ALPHA_DEC, INV_PHI_DEC, SQRT2_MINUS_1_DEC, chi_squared,
                   kronecker, rotation_autocorrelation, serial_correlation,
                   star_discrepancy)

SIZES = (1000, 10000, 100000)
SEED = 2025

generators = {
    "ln C (Cedra)": lambda n: kronecker(n, ALPHA_DEC),
    "1/phi": lambda n: kronecker(n, INV_PHI_DEC),
    "sqrt(2) - 1": lambda n: kronecker(n, SQRT2_MINUS_1_DEC),
    "random (MT19937)": lambda n: (lambda rng: [rng.random() for _ in range(n)])(random.Random(SEED)),
}

print("Star discrepancy D*_N (lower is better) and N * D*_N")
print(f"{'generator':<18}" + "".join(f"{'N=' + str(n):>24}" for n in SIZES))
for name, gen in generators.items():
    row = f"{name:<18}"
    for n in SIZES:
        d = star_discrepancy(gen(n))
        row += f"{d:>14.6f} ({n * d:7.2f})"
    print(row)

print("\nChi-squared, 20 bins, N = 2000 (random reference: ~19)")
for name, gen in generators.items():
    print(f"  {name:<18} {chi_squared(gen(2000)):8.3f}")

print("\nAutocorrelation of {n ln C}: measured (N = 20000) vs exact 1 - 6f(1-f)")
data = kronecker(20000)
for lag in (1, 2, 3, 5, 8, 13, 21, 34, 55, 89):
    print(f"  lag {lag:>3}: {serial_correlation(data, lag):+.3f}   exact {rotation_autocorrelation(lag, ALPHA):+.3f}")
print("High values at Fibonacci lags come from alpha being close to 1/phi.")
print("Consecutive pairs (X_n, X_n+1) lie on two lines: unusable for 2D sampling as is.")
