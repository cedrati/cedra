"""
Diffraction of the Cedra quasicrystal: numerical check of docs/diffraction.md.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Bragg peaks sit at k = {j alpha} (mod 1) with intensity
I(j) = |sin(pi j tau) / (pi j)|^2, the intensities add up to the density tau
(so there is no continuous part), and the autocorrelation is
eta(m) = tau - min(||m alpha||, tau, 1 - tau).
"""

import math

from cedra import (ALPHA, ALPHA_DEC, INV_PHI_DEC, autocorrelation, bragg_intensity,
                   fourier_bohr, kronecker, window_word)

N = 1_000_000
PEAKS = (1, 2, 3, 4, 5, 8, 13, 47)
OFF_LATTICE = (0.3, math.e - 2, math.sqrt(3) - 1)


def analyse(name, window):
    tau = float(window)
    print("\n" + "=" * 72)
    print(f"Window {name}, tau = {tau:.6f}, N = {N:,}")
    print("=" * 72)
    word = window_word(N, window)
    fractions = kronecker(N)
    ok = True

    print(f"{'j':>4} {'k = {j alpha}':>14} {'theory I(j)':>13} {'measured':>13}  check")
    for j in PEAKS:
        theory = bragg_intensity(j, window)
        measured = abs(fourier_bohr(word, fractions, j=j)) ** 2
        # Finite-size leakage from the central peak (k = 0) is about (tau / (pi N ||j alpha||))^2
        distance = min((j * ALPHA) % 1, 1 - (j * ALPHA) % 1)
        leakage = (tau / (math.pi * N * distance)) ** 2
        good = abs(measured - theory) < max(0.02 * theory, 4 * leakage, 1e-9)
        ok &= good
        print(f"{j:>4} {(j * ALPHA) % 1:>14.6f} {theory:>13.4e} {measured:>13.4e}  {'PASS' if good else 'FAIL'}")

    print("\nOff the module Z + alpha Z (should tend to 0 like 1/N^2):")
    for k in OFF_LATTICE:
        print(f"  k = {k:.6f}: |a(k)|^2 = {abs(fourier_bohr(word, fractions, k=k)) ** 2:.2e}")

    total = bragg_intensity(0, window) + 2 * sum(bragg_intensity(j, window) for j in range(1, 200001))
    tail = 1 / (math.pi ** 2 * 200000)
    good = abs(total - tau) < tail
    ok &= good
    print(f"\nTotal Bragg intensity sum_{{|j| <= 2e5}} I(j) = {total:.8f} vs tau = {tau:.8f}"
          f" (tail < {tail:.1e})  {'PASS' if good else 'FAIL'}")

    print("\nAutocorrelation eta(m) = tau - min(||m alpha||, tau, 1 - tau):")
    for m in (1, 2, 3, 13, 47, 1376):
        measured = sum(word[i] and word[i + m] for i in range(N - m)) / (N - m)
        theory = float(autocorrelation(m, window))
        good = abs(measured - theory) < 1e-3
        ok &= good
        print(f"  m = {m:>4}: measured {measured:.5f}, theory {theory:.5f}  {'PASS' if good else 'FAIL'}")
    return ok


if __name__ == "__main__":
    ok = analyse("1/phi (model set)", INV_PHI_DEC)
    ok &= analyse("alpha (Sturmian)", ALPHA_DEC)
    print("\nNote: for the Sturmian word, I(q_k) = sin^2(pi q_k alpha) / (pi q_k)^2 is almost zero:")
    for q in (13, 47, 1376):
        print(f"  I({q}) = {bragg_intensity(q, ALPHA_DEC):.2e}")
    print("\nAll checks passed." if ok else "\nSome checks FAILED.")
