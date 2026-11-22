# Changelog

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
