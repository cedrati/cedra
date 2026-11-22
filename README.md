# Cedra

*An irrational rotation close to the golden ratio, and the temporal quasicrystal it generates*

<br/>

Cedra is the algebraic number

$$C = \sqrt{3} + \sqrt{2} + \sqrt{\tfrac{1}{2}} - 2 = \sqrt{3} + \tfrac{3}{2}\sqrt{2} - 2 \approx 1.853371151128520$$

Its logarithm $\alpha = \ln C \approx 0.617006$ is a transcendental number that happens to sit very close to $1/\varphi \approx 0.618034$. This repository studies the rotation $n \mapsto \{n\alpha\}$ and the one-dimensional quasicrystals (model sets and Sturmian words) it produces. Every number below is reproducible with the scripts of the repository (Python standard library only).

<br/>

## 📐 Algebraic facts (proven)

- **Simplification.** $\sqrt{2} + \sqrt{1/2} = \tfrac{3}{2}\sqrt{2}$, so the definition has three independent terms.
- **Minimal polynomial.** $C$ is algebraic of degree 4: $P(x) = 4x^4 + 32x^3 + 36x^2 - 112x - 167$.
- **Transcendence of $\alpha$.** $C$ is algebraic and $C \neq 0, 1$, so $\ln C$ is transcendental (Hermite, Lindemann, Weierstrass). In particular $\alpha$ is irrational.
- **Continued fraction.** $\alpha = [0; 1, 1, 1, 1, 1, 1, 3, 29, 1, 3, 1, 34, 2, \dots]$, while $1/\varphi = [0; 1, 1, 1, \dots]$.
  Convergent denominators: $1, 1, 2, 3, 5, 8, 13, 47, 1376, 1423, 5645, 7068, \dots$

The continued fraction explains everything that follows: $\alpha$ starts exactly like $1/\varphi$ (Fibonacci denominators up to 13), then the large partial quotients 29 and 34 make it a much better rational approximation target than $1/\varphi$.

> **29/47 is a convergent of $\ln C$.** The approximation $C \approx e^{29/47}$ (relative error $1.5 \times 10^{-5}$) is simply the 7th convergent of the continued fraction. Its accuracy comes from the next partial quotient being large (29), not from a hidden harmony.

<br/>

## 🔮 Temporal quasicrystal

The rotation is coded by a window $W \subset [0, 1)$:

$$S_W = \{\, n \in \mathbb{Z} \mid \{n\alpha\} \in W \,\}$$

where $\{x\} = x - \lfloor x \rfloor$. This is the cut-and-project construction of a 1D model set: lattice $\mathbb{Z}^2$, physical space along $n$, internal space the circle.

| Window | Word | Factor complexity $p(n)$ | Gaps | Density |
|---|---|---|---|---|
| $[0, \alpha)$ | **Sturmian** | $n + 1$ | 1, 2 | $\alpha \approx 0.6170$ |
| $[0, 1/\varphi)$ | rotation coding, **not** Sturmian ($1/\varphi \notin \mathbb{Z} + \alpha\mathbb{Z}$) | $2n$ for $n \ge 3$ | 1, 2 | $1/\varphi \approx 0.6180$ |

Proven properties (general theory of model sets with an interval window):

- **Delone set.** Uniformly discrete and relatively dense; by the three-gap theorem, the gaps take at most three values (here only 1 and 2).
- **Aperiodic** because $\alpha$ is irrational.
- **Pure point diffraction**, with Bragg peaks on $\mathbb{Z} + \alpha\mathbb{Z}$ (regular model set).
- **Density** equal to the window length (Weyl equidistribution).

### Almost periodic, never periodic

Shifting the word by a convergent denominator $q$ changes a bit only when $\{n\alpha\}$ falls within $\|q\alpha\|$ of one of the two window edges, so the proportion of mismatches is exactly $2\|q\alpha\|$:

| $q$ | $2\lVert q\alpha\rVert$ | measured (N = 50 000) |
|---|---|---|
| 47 | $1.41 \times 10^{-3}$ | $1.40 \times 10^{-3}$ |
| 1376 | $1.13 \times 10^{-3}$ | $1.11 \times 10^{-3}$ |
| 1423 | $2.85 \times 10^{-4}$ | $2.88 \times 10^{-4}$ |
| 7068 | $8.0 \times 10^{-6}$ | 0 mismatch on 20 000 bits |

A finite test such as "no period found in the first N bits" is therefore not evidence of aperiodicity: on 20 000 bits the window $1/\varphi$ word looks exactly 7068-periodic. Aperiodicity comes from the irrationality of $\alpha$, which is proven.

### Temporal reading

Reading $n$ as discrete time, $S_W$ is a set of instants that is ordered (deterministic, pure point spectrum) but never repeats. The term "temporal quasicrystal" is also used in physics for driven quantum systems (Floquet time quasicrystals); Cedra is a purely mathematical, discrete-time model with no claim on physical systems.

<br/>

## 📊 Statistics of $X_n = \{n \ln C\}$

### Uniformity and quasi-Monte Carlo

Star discrepancy $D^*_N$ (lower is better), from `cedra_benchmark.py`:

| $\alpha$ | N = 1 000 | N = 10 000 | N = 100 000 |
|---|---|---|---|
| $\ln C$ (Cedra) | 0.00622 | 0.00042 | 0.000148 |
| $1/\varphi$ | **0.00134** | **0.00026** | **0.000030** |
| $\sqrt{2} - 1$ | 0.00169 | 0.00027 | 0.000032 |
| pseudo-random (MT19937) | 0.0362 | 0.0110 | 0.0017 |

Like every Kronecker sequence, Cedra is far more uniform than random sampling. Among Kronecker sequences it is **not** a good choice: because of its large partial quotients, it is about 5 times less uniform than $1/\varphi$ at N = 100 000. The golden ratio is the optimal 1D Kronecker rotation.

About the chi-squared test: a truly random uniform sample gives $\chi^2 \approx 19$ (20 bins). Cedra gives $\chi^2 = 2.44$ at N = 2000. A very low value means *more regular than random*, which is the point of a low-discrepancy sequence, not a sign of "better randomness".

### Autocorrelation

For any rotation, the asymptotic autocorrelation at lag $k$ has the exact value $1 - 6f(1-f)$ with $f = \{k\alpha\}$:

| lag | 1 | 2 | 3 | 5 | 8 | 13 | 21 | 34 | 55 |
|---|---|---|---|---|---|---|---|---|---|
| correlation | -0.418 | -0.076 | +0.239 | +0.533 | +0.641 | +0.876 | +0.754 | +0.872 | +0.637 |

The peaks at Fibonacci lags come from $\alpha$ being close to $1/\varphi$: $\{F_k \alpha\}$ is close to 0 or 1. Consecutive pairs $(X_n, X_{n+1})$ lie on two lines, so the sequence is unusable as is for 2D sampling (use several independent irrationals instead) and must never be used for cryptography.

<br/>

## 🚀 Quick start

```bash
python3 cedra_start.py          # constant, minimal polynomial, continued fraction
python3 cedra_quasicrystal.py   # model set and Sturmian word, all checks computed
python3 cedra_benchmark.py      # discrepancy, chi-squared, autocorrelation
python3 -m unittest test_cedra  # unit tests
```

```python
from cedra import ALPHA_DEC, kronecker, model_set, window_word

kronecker(5)                       # [0.617006..., 0.234012..., ...]
model_set(20)                      # [1, 2, 4, 5, 7, 9, 10, 12, 13, 15, 17, 18, 20]
window_word(20, window=ALPHA_DEC)  # Sturmian word
```

Constants are computed with `decimal` at 60 digits: with plain floats, $\{n\alpha\}$ is off by $7 \times 10^{-8}$ at $n = 10^9$ and by $10^{-5}$ at $n = 10^{12}$.

`cedra_biquadratic.py` is a side utility: a scale-invariant solver for $ax^4 + bx^2 + c = 0$.

<br/>

## 🔬 Status of the claims

| Statement | Status |
|---|---|
| $C$ algebraic of degree 4, minimal polynomial above | ✅ Proven |
| $\ln C$ transcendental | ✅ Proven (Lindemann, Weierstrass) |
| $S_W$ aperiodic Delone set, pure point spectrum | ✅ Proven (model set theory) |
| Window $[0, \alpha)$ gives a Sturmian word | ✅ Proven |
| Mismatch rate $2\lVert q\alpha\rVert$ at near-periods | ✅ Proven, verified numerically |
| Fibonacci-lag correlations | ✅ Exact formula, explained by $\alpha \approx 1/\varphi$ |
| Better than $1/\varphi$ for quasi-Monte Carlo | ❌ False, about 5x worse |
| Physical interpretation (dimensions, Planck scale) | ❌ Removed: numerical coincidences depending on units |

<br/>

## 🧭 Open directions

- Characterise the Sturmian word with window $\alpha$ through its S-adic representation, read directly from the continued fraction $[0; 1^6, 3, 29, 1, 3, 1, 34, \dots]$.
- Study the diffraction intensities: Bragg peaks at $k = m + j\alpha$ ($m, j \in \mathbb{Z}$), with amplitude $\left|\frac{\sin(\pi j \tau)}{\pi j}\right|$ for a window of length $\tau$ (largest peaks at the Fibonacci-like $j$).
- Is the continued fraction of $\ln C$ "typical" (Khinchin, Lévy constants) like almost every real number? Nothing is known in general for logarithms of algebraic numbers.

<br/>

## License

CC BY-NC 4.0, see [LICENSE](LICENSE).
