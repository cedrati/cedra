"""
Littlewood products for the Cedra pair (ln C1, ln|C2|).

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

The Littlewood conjecture states that inf_n n ||n a|| ||n b|| = 0 for every
pair (a, b). It is open in general. This script gives two kinds of evidence
for the Cedra pair a = ln C1, b = ln|C2|:

  1. brute force records of n ||n a|| ||n b|| for n <= N;
  2. explicit witnesses n = q_k (convergent denominators of a): since
     q_k ||q_k a|| < 1 / a_{k+1}, a large partial quotient forces a small
     Littlewood product. Computed with `digits` digits, so every value is
     certified (q_k^2 stays far below 10^digits).

Usage: python3 cedra_littlewood.py [N] [digits]   (defaults 2000000, 6000)
"""

from decimal import Decimal, ROUND_FLOOR, localcontext
import math
import sys

from cedra import LOG_CONJUGATES_DEC

if hasattr(sys, "set_int_max_str_digits"):   # Python 3.11+: allow str() of huge q_k
    sys.set_int_max_str_digits(0)

PLASTIC = Decimal("1.324717957244746025960908854478097340734404056901733364534")


def brute_force_records(a, b, n_max):
    """Successive minima of n ||n a|| ||n b|| (float arithmetic, fine for n <= 1e8)."""
    a, b = float(a), float(b)
    best, records = 1.0, []
    for n in range(1, n_max + 1):
        x, y = n * a, n * b
        p = n * abs(x - round(x)) * abs(y - round(y))
        if p < best:
            best = p
            records.append((n, p))
    return records


def convergent_witnesses(digits):
    """Records of q_k ||q_k a|| ||q_k b|| along the convergents of a = ln C1, b = ln|C2|."""
    with localcontext() as ctx:
        ctx.prec = digits + 100
        s2, s3 = Decimal(2).sqrt(), Decimal(3).sqrt()
        a = (s3 + 3 * s2 / 2 - 2).ln()
        b = abs(s3 - 3 * s2 / 2 - 2).ln()

        def dist(x):
            f = x - x.to_integral_value(rounding=ROUND_FLOOR)
            return min(f, 1 - f)

        x, q_prev, q, k = a, 0, 1, 0
        best, records = None, []
        while True:
            quotient = int(x.to_integral_value(rounding=ROUND_FLOOR))
            if k >= 1:
                q_prev, q = q, quotient * q + q_prev
            if 2 * q.bit_length() * math.log10(2) > digits:
                return records
            if k >= 1:
                product = q * dist(q * a) * dist(q * b)
                if best is None or product < best:
                    best = product
                    n_digits = len(str(q))
                    records.append((k, n_digits, float(product)))
            x = 1 / (x - quotient)
            k += 1


def main():
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    digits = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    lc = LOG_CONJUGATES_DEC

    print(f"Brute force records of n ||n a|| ||n b|| up to n = {n_max:,}")
    pairs = {
        "Cedra (ln C1, ln|C2|)": (lc[0], lc[1]),
        "R2 (plastic number)": (1 / PLASTIC, 1 / PLASTIC ** 2),
        "(sqrt 2, sqrt 3)": (Decimal(2).sqrt(), Decimal(3).sqrt()),
    }
    for name, (a, b) in pairs.items():
        rec = brute_force_records(a, b, n_max)
        n, p = rec[-1]
        print(f"  {name:<24} best {p:.3e} at n = {n:,}")

    print(f"\nExplicit witnesses along the convergents of ln C ({digits} digits)")
    print(f"  {'k':>5} {'digits of n = q_k':>18} {'n ||n a|| ||n b||':>20}")
    for k, n_digits, product in convergent_witnesses(digits):
        print(f"  {k:>5} {n_digits:>18} {product:>20.3e}")
    print("\nEvidence, not a proof: Littlewood holds for this pair if ln C has unbounded")
    print("partial quotients, which is expected but not known.")


if __name__ == "__main__":
    main()
