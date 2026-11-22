"""
Statistics of the continued fraction of ln C compared with a typical real number.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

For almost every real number (Lebesgue), the partial quotients follow the
Gauss-Kuzmin law, their geometric mean tends to Khinchin's constant and
q_n^(1/n) tends to Levy's constant. Nothing is proven for ln C; this script
only measures it.

Usage: python3 cedra_continued_fraction.py [digits]   (default 6000, about 10 s)
"""

from decimal import Decimal, ROUND_FLOOR, localcontext
import math
import sys

KHINCHIN = 2.6854520010653064
LEVY = math.exp(math.pi ** 2 / (12 * math.log(2)))


def continued_fraction_ln_c(digits):
    """Partial quotients of ln C that are certified by `digits` digits of precision."""
    with localcontext() as ctx:
        ctx.prec = digits + 50
        c = Decimal(3).sqrt() + Decimal(3) / 2 * Decimal(2).sqrt() - 2
        x = c.ln()
        quotients, q_prev, q = [], 0, 1
        while True:
            a = int(x.to_integral_value(rounding=ROUND_FLOOR))
            quotients.append(a)
            if len(quotients) > 1:
                q_prev, q = q, a * q + q_prev
            # |alpha - p/q| < 1/q^2: stop while q^2 is still far below 10^digits
            if 2 * q.bit_length() * math.log10(2) > digits - 50:
                return quotients, q
            x = 1 / (x - a)


def main():
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
    cf, q = continued_fraction_ln_c(digits)
    terms = cf[1:]
    n = len(terms)

    print(f"ln C = [0; {', '.join(map(str, terms[:20]))}, ...]")
    print(f"{n} partial quotients from {digits} digits\n")

    geo = math.exp(sum(math.log(a) for a in terms) / n)
    levy = math.exp(q.bit_length() * math.log(2) / n)
    print(f"Geometric mean     {geo:.4f}   Khinchin constant {KHINCHIN:.4f}")
    print(f"q_n^(1/n)          {levy:.4f}   Levy constant     {LEVY:.4f}\n")

    print("Frequency of each partial quotient vs Gauss-Kuzmin law")
    for k in range(1, 8):
        gk = -math.log2(1 - 1 / (k + 1) ** 2)
        print(f"  a = {k}: {terms.count(k) / n:.4f}   expected {gk:.4f}")

    print("\nLargest partial quotients (value, index):")
    print(" ", sorted(((a, i + 1) for i, a in enumerate(terms)), reverse=True)[:5])
    print(f"\nProbability that a typical number starts with six 1s: {(math.log2(4 / 3)) ** 6:.2%} (roughly)")


if __name__ == "__main__":
    main()
