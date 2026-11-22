"""
Cedra temporal quasicrystal: cut-and-project construction.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Two codings of the rotation n -> {n alpha}, alpha = ln C:

  * model set with window 1/phi:  S = {n | {n alpha} < 1/phi}
    aperiodic Delone set, factor complexity 2n for n >= 3 (not Sturmian,
    since 1/phi is not in Z + alpha Z)
  * Sturmian word with window alpha: s(n) = 1 if {n alpha} < alpha
    aperiodic, factor complexity exactly n + 1

Both words are aperiodic because alpha is irrational, yet almost periodic
with near-periods 47, 1376, 1423, 5645, 7068... (convergent denominators
of ln C). A finite "no period found" test is therefore meaningless: on
20000 bits the window-1/phi word looks exactly 7068-periodic.

Every check below is computed; nothing is printed as passed unless it is.
"""

from cedra import (ALPHA, ALPHA_DEC, INV_PHI, INV_PHI_DEC, convergent_denominators,
                   distance_to_integer, factor_complexity, model_set, window_word)

N = 50000
LENGTHS = (5, 10, 50, 100)
NEAR_PERIODS = [q for q in convergent_denominators(ALPHA_DEC, 12) if 40 < q < N // 4]


def report(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


def analyse(title, window, window_float, expected_complexity):
    print("\n" + "=" * 64)
    print(title)
    print("=" * 64)

    points = model_set(N, window)
    word = window_word(N, window)
    print(f"First points: {points[:20]} ...")
    print(f"First bits:   {''.join(map(str, word[:60]))}")

    gaps = sorted({b - a for a, b in zip(points, points[1:])})
    density = len(points) / N
    complexities = {k: factor_complexity(word, k) for k in LENGTHS}

    results = [
        report("Delone set", len(gaps) <= 3,
               f"gaps take the values {gaps} (three-gap theorem allows at most 3)"),
        report("Density", abs(density - window_float) < 5 / N ** 0.5,
               f"{density:.6f} vs window length {window_float:.6f}"),
        report("Factor complexity",
               all(complexities[k] == expected_complexity(k) for k in LENGTHS),
               ", ".join(f"p({k})={v}" for k, v in complexities.items())),
    ]

    # Aperiodicity is a theorem (alpha is irrational), but the word is *almost*
    # periodic: shifting by a convergent denominator q changes a bit only when
    # {n alpha} falls within ||q alpha|| of one of the two window edges.
    print("  Near-periods q (convergent denominators of ln C):")
    for q in NEAR_PERIODS:
        rate = sum(word[i] != word[i + q] for i in range(N - q)) / (N - q)
        predicted = 2 * float(distance_to_integer(q * ALPHA_DEC))
        if predicted * (N - q) < 10:
            print(f"    [skip] q = {q:>5}: too few expected mismatches for N = {N}")
            continue
        ok = abs(rate / predicted - 1) < 0.15
        results.append(ok)
        print(f"    [{'PASS' if ok else 'FAIL'}] q = {q:>5}: mismatch rate {rate:.2e}, "
              f"predicted 2||q alpha|| = {predicted:.2e}")
    return all(results)


if __name__ == "__main__":
    ok_model = analyse("MODEL SET, window 1/phi  (expected complexity 2n)",
                       INV_PHI_DEC, INV_PHI, lambda k: 2 * k)
    ok_sturm = analyse("STURMIAN WORD, window alpha = ln C  (expected complexity n + 1)",
                       ALPHA_DEC, ALPHA, lambda k: k + 1)

    print("\n" + "=" * 64)
    print("All checks passed." if ok_model and ok_sturm else "Some checks FAILED.")
