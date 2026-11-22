"""
Cedra constant: quick start.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Shows the definition of the constant, its algebraic facts, and why the
rotation number alpha = ln C behaves almost like 1/phi.
"""

from decimal import Decimal

from cedra import (ALPHA, ALPHA_DEC, CEDRA, CEDRA_DEC, INV_PHI, INV_PHI_DEC,
                   MINIMAL_POLYNOMIAL, continued_fraction, eval_poly, kronecker)

print(f"C = sqrt(3) + sqrt(2) + sqrt(1/2) - 2 = {CEDRA_DEC:.30f}")

# sqrt(2) + sqrt(1/2) = (3/2) sqrt(2): the definition really has three terms
simplified = Decimal(3).sqrt() + Decimal(3) / 2 * Decimal(2).sqrt() - 2
print(f"C = sqrt(3) + (3/2) sqrt(2) - 2      (difference: {abs(CEDRA_DEC - simplified):.1e})")

residual = eval_poly(MINIMAL_POLYNOMIAL, CEDRA_DEC)
print("\nMinimal polynomial 4x^4 + 32x^3 + 36x^2 - 112x - 167")
print(f"P(C) = {residual:.1e}  -> C is algebraic of degree 4")

print(f"\nalpha = ln C = {ALPHA_DEC:.30f}")
print("C is algebraic and C != 0, 1, so ln C is transcendental (Lindemann-Weierstrass).")
print("In particular alpha is irrational and {n alpha} is equidistributed (Weyl).")

print(f"\n1/phi         = {INV_PHI_DEC:.30f}")
print(f"alpha - 1/phi = {ALPHA - INV_PHI:+.6f}")
print(f"Continued fraction of alpha: {continued_fraction(ALPHA_DEC, 15)}")
print(f"Continued fraction of 1/phi: {continued_fraction(INV_PHI_DEC, 15)}")
print("alpha starts like 1/phi (six 1s), which explains the Fibonacci-lag correlations;")
print("the later large quotients (29, 34) make it less uniform than 1/phi itself.")

print("\nKronecker sequence X_n = {n ln C}, first 10 values:")
for n, x in enumerate(kronecker(10), start=1):
    print(f"  X{n:<2} = {x:.6f}")

print(f"\n(float check: C = {CEDRA}, alpha = {ALPHA})")
