# Diffraction of the Cedra quasicrystal

*Olivier Cedrati, November 2026*

A quasicrystal is recognised by its diffraction pattern: sharp Bragg peaks and nothing else. This note computes the diffraction of the Cedra model sets **exactly**: the position and intensity of every peak, together with an elementary proof that there is no continuous part. The only tools are the arc lemma of [near-periods.md](near-periods.md) and the uniqueness of Fourier coefficients on the circle. Numbers are reproduced by [`cedra_diffraction.py`](../cedra_diffraction.py) and covered by [`test_cedra.py`](../test_cedra.py).

## Setting

$\alpha = \ln C$ is irrational, $0 < \tau < 1$, $W = [0, \tau)$, and

$$\Lambda = \{\, n \in \mathbb{Z} \mid \{n\alpha\} \in W \,\}.$$

The **autocorrelation** of $\Lambda$ is the sequence

$$\eta(m) = \lim_{N\to\infty} \frac{\#\{ |n| \le N : n \in \Lambda,\ n + m \in \Lambda\}}{2N+1}, \qquad m \in \mathbb{Z}.$$

Since $\Lambda \subset \mathbb{Z}$, the autocorrelation measure $\sum_m \eta(m)\,\delta_m$ lives on $\mathbb{Z}$, and its Fourier transform (the **diffraction measure**) is 1-periodic. It is therefore described by a positive measure $\nu$ on the circle $\mathbb{T} = \mathbb{R}/\mathbb{Z}$, the one whose Fourier coefficients are the $\eta(m)$:

$$\eta(m) = \int_{\mathbb{T}} e^{2i\pi k m}\, d\nu(k).$$

$\nu$ exists by Herglotz's theorem, since $\eta$ is positive definite, and it is unique.

## Theorem C (exact diffraction)

$$\boxed{\ \nu = \sum_{j \in \mathbb{Z}} I(j)\, \delta_{\{j\alpha\}}, \qquad I(j) = \frac{\sin^2(\pi j \tau)}{\pi^2 j^2}\ (j \ne 0), \quad I(0) = \tau^2\ }$$

In words:
- the diffraction is **pure point** (Bragg peaks only), supported on $\mathbb{Z} + \alpha\mathbb{Z}$;
- the peak at $k = j\alpha$ (mod 1) has intensity $\sin^2(\pi j\tau)/(\pi j)^2$;
- the total intensity of one period is $\sum_j I(j) = \tau$, the density of $\Lambda$.

### Proof

1. **The autocorrelation is explicit.** By Weyl's theorem, $\eta(m) = \lambda\big(W \cap (W - m\alpha)\big)$. By the arc lemma of the near-periods note, this equals
   $$\eta(m) = f(m\alpha), \qquad f(t) = \lambda\big(W \cap (W - t)\big) = \tau - \min(\lVert t\rVert, \tau, 1-\tau).$$
   (Equivalently: the mismatch density of Theorem A is $2(\tau - \eta(m))$.)

2. **Fourier series of $f$.** $f$ is the correlation $\int 1_W(x) 1_W(x+t)\,dx$, so its Fourier coefficients are $|\widehat{1_W}(j)|^2$, with
   $$\widehat{1_W}(j) = \int_0^\tau e^{-2i\pi jx}\,dx = \frac{1 - e^{-2i\pi j\tau}}{2i\pi j}, \qquad |\widehat{1_W}(j)|^2 = \frac{\sin^2(\pi j\tau)}{\pi^2 j^2} = I(j).$$
   Since $\sum_j I(j) < \infty$, the series converges uniformly: $f(t) = \sum_j I(j)\, e^{2i\pi jt}$ for every $t$.

3. **Identification.** Hence $\eta(m) = f(m\alpha) = \sum_j I(j)\, e^{2i\pi j m\alpha}$. These are the Fourier coefficients of the discrete measure $\mu = \sum_j I(j)\,\delta_{\{j\alpha\}}$; its atoms are distinct because $\alpha$ is irrational. A measure on $\mathbb{T}$ is determined by its Fourier coefficients, so $\nu = \mu$.

4. **Total mass.** $\nu(\mathbb{T}) = \eta(0) = \tau$, and also $\sum_j I(j) = \lVert 1_W\rVert_2^2 = \tau$ by Parseval. $\square$

The same computation by the Fourier–Bohr coefficients gives $a(j\alpha) = \widehat{1_W}(j)$, because $e^{-2i\pi j\alpha n} = e^{-2i\pi j\{n\alpha\}}$. So $\nu(\{j\alpha\}) = |a(j\alpha)|^2$, as predicted by the general theory of model sets (Hof, Schlottmann), which is not needed here.

## Numerical check ($N = 10^6$)

| $j$ | peak $k = \{j\alpha\}$ | $I(j)$, window $1/\varphi$ | measured | $I(j)$, window $\alpha$ | measured |
|---|---|---|---|---|---|
| 1 | 0.617006 | $8.8016 \times 10^{-2}$ | $8.8016 \times 10^{-2}$ | $8.8236 \times 10^{-2}$ | $8.8237 \times 10^{-2}$ |
| 2 | 0.234012 | $1.1558 \times 10^{-2}$ | $1.1558 \times 10^{-2}$ | $1.1395 \times 10^{-2}$ | $1.1395 \times 10^{-2}$ |
| 3 | 0.851019 | $2.2041 \times 10^{-3}$ | $2.2041 \times 10^{-3}$ | $2.2913 \times 10^{-3}$ | $2.2912 \times 10^{-3}$ |
| 4 | 0.468025 | $6.2842 \times 10^{-3}$ | $6.2841 \times 10^{-3}$ | $6.2689 \times 10^{-3}$ | $6.2688 \times 10^{-3}$ |
| 5 | 0.085031 | $3.1662 \times 10^{-4}$ | $3.1665 \times 10^{-4}$ | $2.8240 \times 10^{-4}$ | $2.8242 \times 10^{-4}$ |
| 13 | 0.021081 | $6.99 \times 10^{-6}$ | $6.98 \times 10^{-6}$ | $2.63 \times 10^{-6}$ | $2.62 \times 10^{-6}$ |

- Off the module $\mathbb{Z} + \alpha\mathbb{Z}$ (for example $k = 0.3$, $e - 2$, $\sqrt{3} - 1$), the measured intensity is below $5 \times 10^{-11}$.
- The autocorrelation formula matches at every lag tested ($m = 1, 2, 3, 13, 47, 1376$) to 5 digits.
- Finite-size effect: the peak $j = 47$ sits at distance $\lVert 47\alpha\rVert = 7 \times 10^{-4}$ from the huge central peak ($I(0) = \tau^2 \approx 0.38$). At size $N$, the central peak leaks about $(\tau / \pi N \lVert j\alpha\rVert)^2$ onto it. Resolving such a peak needs $N \gg 1/\lVert q_k\alpha\rVert$, so the convergent denominators appear again.

## Remarks

**Peaks and near-periods are dual.** The near-period $q_k$ of the near-periods note shows up here as the peaks $\{q_k\alpha\}$, $\{2q_k\alpha\}, \dots$ crowding at distance $\lVert q_k\alpha\rVert$ from $k = 0$. A sharp near-period (large next quotient, like 29 after 47) corresponds to peaks that are hard to separate from the origin.

**Extinctions of the Sturmian word.** For $\tau = \alpha$, $I(j) = \sin^2(\pi j\alpha)/(\pi j)^2$ vanishes as $j\alpha$ approaches an integer. At the convergents, $I(q_k) \approx \lVert q_k\alpha\rVert^2 / q_k^2$:

| $q_k$ | 13 | 47 | 1376 |
|---|---|---|---|
| $I(q_k)$, Sturmian | $2.6 \times 10^{-6}$ | $2.3 \times 10^{-10}$ | $1.7 \times 10^{-13}$ |

The Sturmian word thus has almost-extinct peaks exactly where the model set with window $1/\varphi$ still has visible ones.

**The brightest peaks.** $I(j) \le 1/(\pi j)^2$, so the pattern is dominated by $j = 1$ at $k \approx 0.617$, then $j = 2$ at $k \approx 0.234$, and $j = 4$ at $k \approx 0.468$. The peak $j = 3$ is dimmed because $3\tau$ is close to an integer.
