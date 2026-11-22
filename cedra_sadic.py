"""
S-adic structure of the Cedra Sturmian word: numerical check of docs/sturmian-sadic.md.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

The Sturmian word w(n) = 1 iff {n ln C} < ln C (n >= 1) equals 0 c_alpha, where

    c_alpha = sigma_1^5 o sigma_3 o sigma_29 o sigma_1 o sigma_3 o sigma_1 o sigma_34 o ... (1)
    sigma_a : 1 -> 1^a 0,  0 -> 1         (sigma_1 = Fibonacci substitution)
"""

from cedra import (ALPHA_DEC, INV_PHI_DEC, continued_fraction, convergent_denominators,
                   sadic_word, window_word)

N = 1_000_000


def standard_words(alpha, count):
    """s_{-1} = 0, s_0 = 1, s_k = s_{k-1}^{a_{k+1}} s_{k-2} for the quotients a_2, a_3, ..."""
    cf = continued_fraction(alpha, count + 2)
    words = ['0', '1']
    for a in cf[2:]:
        words.append(words[-1] * a + words[-2])
    return words[1:]


def common_prefix(u, v):
    for i, (x, y) in enumerate(zip(u, v)):
        if x != y:
            return i
    return min(len(u), len(v))


if __name__ == "__main__":
    ok = True
    cf = continued_fraction(ALPHA_DEC, 14)
    print(f"ln C = [0; {', '.join(map(str, cf[1:]))}, ...]")
    print("c_alpha = " + " o ".join(f"sigma_{a}" for a in cf[2:]) + " o ... (1)")

    w = ''.join(map(str, window_word(N + 1, ALPHA_DEC)))
    c = sadic_word(N)
    good = w[0] == '0' and w[1:] == c
    ok &= good
    print(f"\n[{'PASS' if good else 'FAIL'}] w = 0 c_alpha on the first {N:,} letters")

    print("\nStandard words s_k (prefixes of c_alpha), lengths = convergent denominators:")
    words = standard_words(ALPHA_DEC, 10)
    dens = convergent_denominators(ALPHA_DEC, 12)[1:]
    for k, s in enumerate(words):
        is_prefix = c.startswith(s[:-2])
        good = len(s) == dens[k] and is_prefix
        ok &= good
        shown = s if len(s) <= 40 else s[:37] + "..."
        print(f"  [{'PASS' if good else 'FAIL'}] |s_{k}| = {len(s):>7}  {shown}")

    fib = sadic_word(N, INV_PHI_DEC)
    shared = common_prefix(c, fib)
    good = shared == 2 * 13 + 8 - 2
    ok &= good
    print(f"\n[{'PASS' if good else 'FAIL'}] common prefix with the Fibonacci word: {shared} letters "
          f"(predicted 2*13 + 8 - 2 = 32)")
    print(f"  Cedra    : {c[:40]}")
    print(f"  Fibonacci: {fib[:40]}")
    print(f"             {' ' * shared}^")
    print("\nAll checks passed." if ok else "\nSome checks FAILED.")
