"""
Near-periods of the Cedra quasicrystal: numerical check of docs/near-periods.md.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Theorem A  density of {n : w(n) != w(n + q)} = 2 min(||q alpha||, tau, 1 - tau)
Theorem B  Sturmian word (window alpha), q = q_k convergent denominator:
           if q_k alpha - p_k > 0 the first mismatch is at n = q_{k+1},
           otherwise it is at n = 1.
"""

from decimal import Decimal

from cedra import (ALPHA_DEC, INV_PHI_DEC, continued_fraction, convergent_denominators,
                   distance_to_integer, mismatch_density, mismatch_positions, window_word)

N = 1_000_000
CONVERGENTS = convergent_denominators(ALPHA_DEC, 14)   # 1, 1, 2, 3, 5, ..., 498982
SHIFTS = sorted(set(CONVERGENTS[2:]) | {4, 10, 100, 1000})


def theorem_a(name, window):
    print(f"\nTheorem A, window {name} (N = {N:,})")
    print(f"{'q':>8} {'||q alpha||':>12} {'predicted':>11} {'measured':>11} {'mismatches':>11}  check")
    word = window_word(N, window)
    all_ok = True
    for q in SHIFTS:
        predicted = float(mismatch_density(q, window))
        count = len(mismatch_positions(word, q))
        measured = count / (N - q)
        # Weyl error is O(D*_N): only compare when enough mismatches are expected
        if predicted * (N - q) < 30:
            status = "skip (too few expected)"
        else:
            ok = abs(measured / predicted - 1) < 0.05
            all_ok &= ok
            status = "PASS" if ok else "FAIL"
        print(f"{q:>8} {float(distance_to_integer(q * ALPHA_DEC)):>12.3e} "
              f"{predicted:>11.4e} {measured:>11.4e} {count:>11}  {status}")
    return all_ok


def theorem_b():
    print(f"\nTheorem B, Sturmian word (window alpha), N = {N:,}")
    print(f"{'q_k':>8} {'sign of q_k alpha - p_k':>24} {'predicted first':>16} {'measured first':>15}  check")
    word = window_word(N, ALPHA_DEC)
    alpha = Decimal(ALPHA_DEC)
    all_ok = True
    qs = CONVERGENTS[2:]
    for q, q_next in zip(qs, qs[1:]):
        above = (q * alpha) % 1 < Decimal("0.5")      # q alpha slightly above an integer
        predicted = q_next if above else 1
        if predicted + q > N:
            continue
        positions = mismatch_positions(word[:predicted + q + 1], q)
        first = positions[0] if positions else None
        ok = first == predicted
        all_ok &= ok
        print(f"{q:>8} {'+' if above else '-':>24} {predicted:>16} {str(first):>15}  {'PASS' if ok else 'FAIL'}")
    return all_ok


def sharpness():
    print("\nSharpness q_k ||q_k alpha|| = 1 / (alpha_{k+1} + q_{k-1} / q_k)")
    cf = continued_fraction(ALPHA_DEC, 16)
    qs = convergent_denominators(ALPHA_DEC, 16)
    for k in range(2, 13):
        product = qs[k] * distance_to_integer(qs[k] * ALPHA_DEC)
        print(f"  q_{k:<2} = {qs[k]:>7}   next quotient a_{k+1} = {cf[k + 1]:>3}   "
              f"q ||q alpha|| = {float(product):.4f}")
    print("  Large partial quotients (29, 34) give very sharp near-periods (47, 7068),")
    print("  and the same quotients are what make the discrepancy of {n ln C} worse than 1/phi.")


if __name__ == "__main__":
    ok = theorem_a("1/phi", INV_PHI_DEC)
    ok &= theorem_a("alpha (Sturmian)", ALPHA_DEC)
    ok &= theorem_b()
    sharpness()
    print("\nAll checks passed." if ok else "\nSome checks FAILED.")
