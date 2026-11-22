# Cedra

**An irrational rotation close to the golden ratio, and the one-dimensional quasicrystals it generates.**

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](LICENSE)
![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)
![Dependencies: none](https://img.shields.io/badge/dependencies-none-brightgreen.svg)

Cedra is the algebraic number

$$C = \sqrt{3} + \sqrt{2} + \sqrt{\tfrac{1}{2}} - 2 = \sqrt{3} + \tfrac{3}{2}\sqrt{2} - 2 \approx 1.853371151128520.$$

Its logarithm $\alpha = \ln C \approx 0.617006$ is transcendental and lies very close to $1/\varphi \approx 0.618034$. This repository studies the rotation $n \mapsto \{n\alpha\}$ and the structures it generates: model sets, Sturmian words, their diffraction, and their extension to several dimensions through the Galois conjugates of $C$.

Every result is either **proven** in the notes or **measured** by a script, and all numbers are reproducible with the Python standard library alone.

---

## Contents

- [Main results](#main-results)
- [Background](#background)
- [Repository layout](#repository-layout)
- [Usage](#usage)
- [Numerical highlights](#numerical-highlights)
- [Scope and limitations](#scope-and-limitations)
- [Open questions](#open-questions)
- [Citation](#citation)
- [License](#license)

---

## Main results

| # | Statement | Status | Reference |
|---|---|---|---|
| 1 | $C$ is algebraic of degree 4, with minimal polynomial $4x^4 + 32x^3 + 36x^2 - 112x - 167$ | Proven | [`cedra_start.py`](cedra_start.py) |
| 2 | $\alpha = \ln C$ is transcendental (Lindemann, Weierstrass) | Proven | [Background](#background) |
| 3 | **Theorem A.** For every shift $q$, the density of $n$ with $w(n) \neq w(n+q)$ is $2\min(\lVert q\alpha\rVert, \tau, 1-\tau)$ | Proven, verified to $N = 10^6$ | [near-periods](docs/near-periods.md) |
| 4 | **Theorem B.** In the Sturmian word, the near-period $q_k$ fails first at $n = q_{k+1}$ or at $n = 1$, depending on the sign of $q_k\alpha - p_k$ | Proven, verified | [near-periods](docs/near-periods.md) |
| 5 | **Theorem C.** The diffraction is pure point: $\nu = \sum_{j} \frac{\sin^2(\pi j\tau)}{\pi^2 j^2}\,\delta_{\{j\alpha\}}$, with an elementary proof | Proven, verified | [diffraction](docs/diffraction.md) |
| 6 | The Sturmian word is $\sigma_1^5 \circ \sigma_3 \circ \sigma_{29} \circ \sigma_1 \circ \sigma_3 \circ \cdots (1)$, where $\sigma_1$ is the Fibonacci substitution; it agrees with the Fibonacci word on exactly 32 letters | Proven, verified | [S-adic](docs/sturmian-sadic.md) |
| 7 | The Sturmian word is neither a substitution fixed point nor $k$-automatic | Proven (from 2) | [S-adic](docs/sturmian-sadic.md) |
| 8 | **Theorem D.** $1, \ln\lvert C_1\rvert, \dots, \ln\lvert C_4\rvert$ (Galois conjugates) are linearly independent over $\mathbb{Q}$, via the totally split prime 167 | Proven | [higher dimensions](docs/higher-dimensions.md) |
| 9 | An explicit integer with 2 117 digits gives the Littlewood product $n\lVert n\ln C_1\rVert\lVert n\ln\lvert C_2\rvert\rVert < 1.5 \times 10^{-7}$ | Certified computation | [higher dimensions](docs/higher-dimensions.md#littlewood-products) |
| 10 | $\{n \ln C\}$ is about 5 times *less* uniform than $\{n/\varphi\}$ at $N = 10^5$ | Measured | [`cedra_benchmark.py`](cedra_benchmark.py) |
| 11 | The continued fraction of $\ln C$ looks statistically typical (Khinchin, Lévy, Gauss–Kuzmin) | Measured, not provable with current methods | [near-periods](docs/near-periods.md) |

## Background

**Continued fraction.** Everything is governed by

$$\alpha = [0;\ 1, 1, 1, 1, 1, 1,\ 3,\ 29,\ 1, 3, 1,\ 34,\ 2, \dots]$$

with convergent denominators $1, 1, 2, 3, 5, 8, 13, 47, 1376, 1423, 5645, 7068, \dots$

- The opening run of six 1s is shared with $1/\varphi = [0; 1, 1, 1, \dots]$. It explains the "golden ratio flavour": Fibonacci-lag correlations, and a Sturmian word that starts like the Fibonacci word.
- The large partial quotients 29 and 34 explain three facts at once:
  - the approximation $C \approx e^{29/47}$, since $29/47$ is the 7th convergent;
  - unusually sharp near-periods at 47 and 7068;
  - a discrepancy worse than that of $1/\varphi$.

**Transcendence.** $C$ is algebraic and $C \neq 0, 1$, so $\ln C$ is transcendental by the Lindemann–Weierstrass theorem. In particular $\alpha$ is irrational, and $(\{n\alpha\})$ is equidistributed (Weyl).

**Quasicrystal.** For a window $W = [0, \tau)$,

$$S_W = \{\, n \in \mathbb{Z} \mid \{n\alpha\} \in W \,\}$$

is a one-dimensional cut-and-project (model) set: an aperiodic Delone set with gaps 1 and 2 and density $\tau$. Reading $n$ as discrete time gives a *temporal quasicrystal*: deterministic and ordered, but never periodic.

| Window | Word | Factor complexity | Density |
|---|---|---|---|
| $[0, \alpha)$ | Sturmian | $n + 1$ | $\alpha \approx 0.6170$ |
| $[0, 1/\varphi)$ | rotation coding, not Sturmian | $2n$ for $n \ge 3$ | $1/\varphi \approx 0.6180$ |

## Repository layout

```
cedra.py                     Core library: constants (60-digit decimal), sequences, words, measures
cedra_start.py               Constant, minimal polynomial, continued fraction
cedra_quasicrystal.py        Model set and Sturmian word, computed checks
cedra_near_periods.py        Theorems A and B, N = 10^6
cedra_diffraction.py         Theorem C: Bragg peaks, total intensity, autocorrelation
cedra_sadic.py               S-adic construction, comparison with the Fibonacci word
cedra_multidim.py            Theorem D: Galois conjugates, 2D and 4D discrepancy
cedra_littlewood.py          Littlewood products, explicit witnesses along the convergents
cedra_continued_fraction.py  Continued fraction statistics of ln C
cedra_benchmark.py           1D discrepancy, chi-squared and autocorrelation benchmark
cedra_biquadratic.py         Side utility: scale-invariant solver for ax^4 + bx^2 + c = 0
test_cedra.py                Unit tests (36)
docs/                        Mathematical notes with proofs
```

## Usage

Requirements: Python 3.8 or later, no third-party package.

```bash
git clone https://github.com/cedrati/cedra.git
cd cedra
python3 -m unittest test_cedra      # run the test suite
python3 cedra_near_periods.py       # reproduce any note, e.g. Theorems A and B
```

```python
from cedra import ALPHA_DEC, kronecker, model_set, window_word, sadic_word

kronecker(3)                        # [0.617006..., 0.234012..., 0.851019...]
model_set(20)                       # [1, 2, 4, 5, 7, 9, 10, 12, 13, 15, 17, 18, 20]
window_word(10, window=ALPHA_DEC)   # Sturmian word: [0, 1, 0, 1, 1, 0, 1, 0, 1, 1]
sadic_word(13)                      # '1011010110110'
```

Constants are computed with `decimal` at 60 significant digits. With double-precision floats, $\{n\alpha\}$ would be off by $7 \times 10^{-8}$ at $n = 10^9$ and by $10^{-5}$ at $n = 10^{12}$.

## Numerical highlights

**Near-periods** ($N = 10^6$, window $1/\varphi$). The measured mismatch density matches Theorem A:

| shift $q$ | 47 | 1376 | 1423 | 5645 |
|---|---|---|---|---|
| theory $2\lVert q\alpha\rVert$ | $1.415 \times 10^{-3}$ | $1.130 \times 10^{-3}$ | $2.846 \times 10^{-4}$ | $2.765 \times 10^{-4}$ |
| measured | $1.414 \times 10^{-3}$ | $1.130 \times 10^{-3}$ | $2.844 \times 10^{-4}$ | $2.756 \times 10^{-4}$ |

As a consequence, a finite test such as "no period in the first $N$ bits" proves nothing: 20 000 bits of this word are exactly 7068-periodic.

**One-dimensional uniformity**, star discrepancy $D^*_N$ (lower is better):

| sequence | $N = 10^3$ | $N = 10^4$ | $N = 10^5$ |
|---|---|---|---|
| $\{n \ln C\}$ | 0.00622 | 0.00042 | 0.000148 |
| $\{n/\varphi\}$ | **0.00134** | **0.00026** | **0.000030** |
| $\{n(\sqrt{2}-1)\}$ | 0.00169 | 0.00027 | 0.000032 |
| pseudo-random (MT19937) | 0.0362 | 0.0110 | 0.0017 |

**Autocorrelation** of $\{n\alpha\}$, exact value $1 - 6f(1-f)$ with $f = \{k\alpha\}$:

| lag $k$ | 1 | 2 | 3 | 5 | 8 | 13 | 21 | 34 |
|---|---|---|---|---|---|---|---|---|
| correlation | -0.418 | -0.076 | +0.239 | +0.533 | +0.641 | +0.876 | +0.754 | +0.872 |

## Scope and limitations

- **Not a random number generator.** $\{n\alpha\}$ is deterministic, strongly correlated, and its consecutive pairs lie on two lines. It must never be used for cryptography. Its very low chi-squared ($\chi^2 = 2.44$ against about 19 for random data, 20 bins) means *more regular than random*, not "more random".
- **Not the best low-discrepancy choice.** In 1D the golden ratio is optimal and Cedra is about 5 times worse. In 2D and 4D, the Cedra Galois orbit beats random sampling but not R2 or Halton (see [higher dimensions](docs/higher-dimensions.md)).
- **No physical claim.** "Temporal quasicrystal" refers here to a discrete-time mathematical model. The term is also used in physics for driven quantum systems, with which this work makes no connection.
- **Context.** The results on near-periods, diffraction and S-adic structure belong to the classical theory of rotation words and model sets. The notes give self-contained proofs and explicit values for this example; they do not claim new general theorems.

## Open questions

- **Continued fraction of $\ln C$.** It looks typical over 11 393 partial quotients (geometric mean 2.708 against Khinchin's 2.685; $q_n^{1/n}$ = 3.295 against Lévy's 3.276), with a giant quotient $a_{4091} = 73\,852$. Nothing is known in general for logarithms of algebraic numbers.
- **Simultaneous approximation.** Is the pair $(\ln C_1, \ln\lvert C_2\rvert)$ badly approximable? Up to $n = 2 \times 10^6$ the normalised records stay above 0.076. Almost every pair is not badly approximable (Khintchine), so a negative answer is expected, but no proof is known.
- **Littlewood conjecture for the Cedra pair.** It would follow from unbounded partial quotients of $\ln C$. Explicit witnesses go down to $1.5 \times 10^{-7}$ (an integer with 2 117 digits), but a proof is out of reach.

## Citation

```bibtex
@misc{cedrati2025cedra,
  author       = {Cedrati, Olivier},
  title        = {Cedra: an irrational rotation close to the golden ratio and its quasicrystals},
  year         = {2025},
  howpublished = {\url{https://github.com/cedrati/cedra}},
  note         = {Notes and code, revised November 2026}
}
```

## License

[CC BY-NC 4.0](LICENSE). © 2025 Olivier Cedrati.
