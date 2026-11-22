# Changelog

## 2026-11-22: higher dimensions and README overhaul

### Added
- `docs/higher-dimensions.md`: Theorem D, linear independence of 1, ln|C_1|, ..., ln|C_4| over Q (Galois conjugates of C), proved with the totally split prime 167 and Lindemann's theorem; elementary proof in dimension 2; 2D and 4D discrepancy and simultaneous approximation measurements.
- `cedra_multidim.py`: number theory checks and multidimensional benchmarks.
- `cedra.py`: `CONJUGATES_DEC`, `LOG_CONJUGATES_DEC`, `kronecker_nd`, `halton`, `l2_star_discrepancy`.
- 6 new tests (34 total).

### Changed
- README rewritten: summary table of results with status and references, background, repository layout, usage, numerical highlights, scope and limitations, open questions, citation.

## 2026-11-22: diffraction and S-adic notes

### Added
- `docs/diffraction.md`: Theorem C, exact diffraction measure sum_j sin^2(pi j tau)/(pi j)^2 delta_{j alpha} with an elementary proof of pure point spectrum (autocorrelation, Fourier series, uniqueness on the circle).
- `docs/sturmian-sadic.md`: S-adic form sigma_1^5 o sigma_3 o sigma_29 o ... (1), exact 32-letter agreement with the Fibonacci word, proofs that the word is neither a substitution fixed point nor automatic.
- `cedra_diffraction.py`, `cedra_sadic.py`: numerical checks up to N = 10^6.
- `cedra.py`: `autocorrelation`, `bragg_intensity`, `fourier_bohr`, `substitution`, `sadic_word`.
- 7 new tests (28 total).

## 2026-11-22: near-periods note

### Added
- `docs/near-periods.md`: Theorem A (mismatch density 2 min(||q alpha||, tau, 1 - tau) for every shift q) and Theorem B (first failure of a Sturmian near-period at q_{k+1} or 1), with proofs and tables.
- `cedra_near_periods.py`: verification of both theorems up to N = 10^6.
- `cedra_continued_fraction.py`: Khinchin, Levy and Gauss-Kuzmin statistics of ln C; giant partial quotient a_4091 = 73852.
- `cedra.py`: `mismatch_density`, `mismatch_positions`.
- Tests for theorems A and B (21 tests).

## 2026-11-22: mathematical review

### Added
- `cedra.py`: core module, constants at 60 digits with `decimal` (no float drift for large n), discrepancy, complexity, continued fraction, convergents.
- `cedra_benchmark.py`: reproducible comparison of ln C, 1/phi, sqrt(2) - 1 and a pseudo-random generator.
- `test_cedra.py`: 19 unit tests.
- README: continued fraction of ln C, transcendence of ln C, 29/47 as a convergent, near-periods with mismatch rate 2||q alpha||, open directions.

### Fixed
- C * delta = phi and tau_C = 1/phi were identities true for any C (2 sqrt(3) + 3 sqrt(2) - 4 = 2C): removed from the claims.
- The window 1/phi word is not Sturmian (complexity 2n); the Sturmian word uses the window alpha.
- Quasi-Monte Carlo claims: ln C is about 5 times less uniform than 1/phi. Low chi-squared is no longer presented as "better randomness".
- `cedra_quasicrystal.py`: checks are computed instead of hard-coded, and the finite periodicity test (which finds false periods 1423, 7068) is replaced by the near-period analysis.
- `cedra_biquadratic.py`: scale-invariant tolerances (x^4 - 5x^2 + 4 scaled by 1e-20 returned no root), special case counter that always reached 100%, wrong header.
- `images/README.md`: file names now match the actual files.

### Removed
- Planck frequency coincidence, digit-sum numerology, "dimensional signature", unsourced benchmark table, reference to a paper that does not exist.
