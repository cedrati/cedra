# Near-periods of the Cedra quasicrystal

*Olivier Cedrati, November 2026*

The words produced by the Cedra rotation are aperiodic, yet they repeat themselves almost exactly at very specific shifts. This note gives an exact formula for how often a shift fails, explains which shifts are "almost periods", and locates the first failure in the Sturmian case. These are classical facts of the theory of rotation words, stated and proved here in a self-contained way for this example. Every number is reproduced by [`cedra_near_periods.py`](../cedra_near_periods.py) and [`cedra_continued_fraction.py`](../cedra_continued_fraction.py), and covered by [`test_cedra.py`](../test_cedra.py).

## Setting

- $C = \sqrt{3} + \tfrac{3}{2}\sqrt{2} - 2$ and $\alpha = \ln C \approx 0.617006$, which is irrational (indeed transcendental).
- $\{x\} = x - \lfloor x \rfloor$ is the fractional part and $\lVert x \rVert$ the distance from $x$ to the nearest integer.
- For a window length $0 < \tau < 1$, the word $w : \mathbb{N} \to \{0, 1\}$ is
  $$w(n) = 1 \iff \{n\alpha\} \in A = [0, \tau).$$
  With $\tau = 1/\varphi$ this is the model set of the README; with $\tau = \alpha$ it is the Sturmian word.
- $\alpha = [0; 1, 1, 1, 1, 1, 1, 3, 29, 1, 3, 1, 34, 2, \dots]$ with convergents $p_k / q_k$ and denominators
  $q_k = 1, 1, 2, 3, 5, 8, 13, 47, 1376, 1423, 5645, 7068, 245957, 498982, \dots$

## Theorem A (mismatch density)

For every integer $q \geq 1$, the set of $n$ such that $w(n) \neq w(n + q)$ has density

$$\boxed{\ \lim_{N\to\infty} \frac{\#\{1 \le n \le N : w(n) \neq w(n+q)\}}{N} = 2 \min\big(\lVert q\alpha \rVert,\ \tau,\ 1 - \tau\big)\ }$$

**Lemma.** On the circle $\mathbb{T} = \mathbb{R}/\mathbb{Z}$, let $A$ be an arc of length $\tau$ and $t \in \mathbb{T}$ with $\lVert t \rVert = \varepsilon \le 1/2$. Then $\lambda\big(A \,\triangle\, (A + t)\big) = 2\min(\varepsilon, \tau, 1 - \tau)$.

*Proof.* $\lambda(A \triangle (A+t)) = 2\big(\tau - \lambda(A \cap (A+t))\big)$, and two arcs of length $\tau$ whose starting points are $\varepsilon$ apart overlap on $\max(\tau - \varepsilon, 0) + \max(\tau + \varepsilon - 1, 0)$. If $\varepsilon \le \min(\tau, 1-\tau)$ this gives $2\varepsilon$; if $\tau \le \varepsilon \le 1 - \tau$ it gives $2\tau$; if $1 - \tau \le \varepsilon \le \tau$ it gives $2(1-\tau)$. Since $\varepsilon \le 1/2$, these are the three cases of $2\min(\varepsilon, \tau, 1-\tau)$. $\square$

*Proof of Theorem A.* Put $x_n = \{n\alpha\}$ and $\delta = \{q\alpha\}$. Then $\{(n+q)\alpha\} = x_n + \delta$ on the circle, so
$$w(n) \neq w(n+q) \iff x_n \in A \,\triangle\, (A - \delta).$$
This set is a finite union of intervals of measure $2\min(\lVert q\alpha\rVert, \tau, 1-\tau)$ by the lemma. By Weyl's equidistribution theorem ($\alpha$ irrational), the proportion of $n \le N$ with $x_n$ in a finite union of intervals tends to its measure. $\square$

The convergence is quantitative: the set is a union of at most four intervals of $[0,1)$, so the error at rank $N$ is at most $8 D^*_N$, where $D^*_N$ is the star discrepancy of $(x_n)$.

**Corollary (aperiodicity, quantified).** Since $q\alpha \notin \mathbb{Z}$, every shift fails on a set of positive density: $w$ is not periodic, not even eventually periodic. The "best" near-periods are the $q$ that minimise $\lVert q\alpha\rVert$, that is the convergent denominators $q_k$ (best approximation property). They fail with density $2\lVert q_k \alpha\rVert < 2/q_{k+1}$.

### Numerical check ($N = 10^6$)

| $q$ | $2\lVert q\alpha\rVert$ (theory) | window $1/\varphi$ | window $\alpha$ |
|---|---|---|---|
| 13 | $4.216 \times 10^{-2}$ | $4.216 \times 10^{-2}$ | $4.216 \times 10^{-2}$ |
| 47 | $1.415 \times 10^{-3}$ | $1.414 \times 10^{-3}$ | $1.415 \times 10^{-3}$ |
| 1376 | $1.130 \times 10^{-3}$ | $1.130 \times 10^{-3}$ | $1.130 \times 10^{-3}$ |
| 1423 | $2.846 \times 10^{-4}$ | $2.844 \times 10^{-4}$ | $2.854 \times 10^{-4}$ |
| 5645 | $2.765 \times 10^{-4}$ | $2.756 \times 10^{-4}$ | $2.756 \times 10^{-4}$ |
| 7068 | $8.05 \times 10^{-6}$ | 8 mismatches | 9 mismatches |

The clipping by $\min(\tau, 1-\tau)$ is also visible: for $q = 4$, $\lVert 4\alpha\rVert = 0.468 > 1 - 1/\varphi$, and the measured density is $0.76393 = 2(1 - 1/\varphi)$.

## Theorem B (first mismatch of the Sturmian word)

Take $\tau = \alpha$. Let $k \ge 1$ with $q_{k+1} > q_k$.

1. If $q_k\alpha - p_k > 0$, then $w(n) = w(n + q_k)$ for all $1 \le n < q_{k+1}$, and $w(q_{k+1}) \neq w(q_{k+1} + q_k)$.
   In other words, the prefix of length $q_{k+1} + q_k - 1$ is exactly $q_k$-periodic, and not one bit longer.
2. If $q_k\alpha - p_k < 0$, then $w(1) \neq w(1 + q_k)$: the shift fails immediately.

Since the sign of $q_k\alpha - p_k$ alternates with $k$, the shifts $q_k$ alternately hold for a long time and fail at once.

*Proof.* Let $\varepsilon = \lVert q_k\alpha\rVert$. We use the best approximation property: for $1 \le m < q_{k+1}$, $\lVert m\alpha\rVert \ge \varepsilon$, with equality only for $m = q_k$. We also use that $q_{k+1}\alpha - p_{k+1}$ has the sign opposite to $q_k\alpha - p_k$ and that $\lVert q_{k+1}\alpha\rVert < \varepsilon$.

*Case 1*, $\delta = \{q_k\alpha\} = \varepsilon$. A direct check shows that $x \in [0,\alpha)$ and $x + \varepsilon \in [0, \alpha)$ (mod 1) disagree exactly when $x \in [\alpha - \varepsilon, \alpha) \cup [1 - \varepsilon, 1)$.
- $x_n \in [1-\varepsilon, 1)$ means that $n\alpha$ lies *below* an integer at distance at most $\varepsilon$. For $n < q_{k+1}$ this forces $n = q_k$, but $q_k\alpha$ lies above an integer. So there is no such $n$.
- $x_n \in [\alpha - \varepsilon, \alpha)$ means $x_{n-1} \in [1 - \varepsilon, 1)$ for $n \ge 2$ (and $x_1 = \alpha$ is not in it). The same argument excludes every $n - 1 < q_{k+1}$.

So no mismatch occurs for $n < q_{k+1}$. At $n = q_{k+1}$, $q_{k+1}\alpha$ lies below an integer at distance $< \varepsilon$, so $x_{q_{k+1}} \in (1-\varepsilon, 1)$: mismatch.

*Case 2*, $\delta = 1 - \varepsilon$. The disagreement set is $[0, \varepsilon) \cup [\alpha, \alpha + \varepsilon)$, which contains $x_1 = \alpha$. $\square$

### Numerical check

| $q_k$ | sign of $q_k\alpha - p_k$ | predicted first mismatch | measured |
|---|---|---|---|
| 13 | + | 47 | 47 |
| 47 | − | 1 | 1 |
| 1376 | + | 1423 | 1423 |
| 1423 | − | 1 | 1 |
| 5645 | + | 7068 | 7068 |
| 7068 | − | 1 | 1 |
| 245957 | + | 498982 | 498982 |

For the window $1/\varphi$ there is no such clean rule, because the window edge $1/\varphi$ is unrelated to $\alpha$. There the first mismatches are 1330, 1377, 2753, 4176 and 179453 for $q = 47, 1376, 1423, 5645$ and 7068. This is why a "no period in the first N bits" test can be fooled: 20 000 bits of that word are exactly 7068-periodic.

## How sharp is a near-period?

A standard identity of continued fractions gives
$$q_k \lVert q_k\alpha\rVert = \frac{1}{\alpha_{k+1} + q_{k-1}/q_k}, \qquad \alpha_{k+1} = [a_{k+1}; a_{k+2}, \dots].$$
So a near-period $q_k$ is sharp exactly when the *next* partial quotient is large:

| $q_k$ | 13 | **47** | 1376 | 1423 | 5645 | **7068** | 245957 |
|---|---|---|---|---|---|---|---|
| $a_{k+1}$ | 3 | **29** | 1 | 3 | 1 | **34** | 2 |
| $q_k\lVert q_k\alpha\rVert$ | 0.274 | **0.033** | 0.778 | 0.203 | 0.781 | **0.028** | 0.368 |

This explains three facts of the README with a single cause:

- $e^{29/47} \approx C$ is accurate because $29/47$ is followed by the quotient 29;
- the near-periods 47 and 7068 are unusually sharp;
- $\{n \ln C\}$ has a worse discrepancy than $\{n/\varphi\}$, whose quotients are all 1 (the bounded quotients of $1/\varphi$ are exactly what make it optimal).

## Is ln C a typical number?

For almost every real number, the continued fraction follows the Gauss–Kuzmin law, the geometric mean of the quotients tends to Khinchin's constant, and $q_n^{1/n}$ tends to Lévy's constant. Nothing of this is proven for any explicit logarithm, so we can only measure it. With 11 393 quotients (12 000 digits; 5 764 with the default 6 000 digits of the script):

| | ln C | typical value |
|---|---|---|
| geometric mean of $a_n$ | 2.708 | 2.685 (Khinchin) |
| $q_n^{1/n}$ | 3.295 | 3.276 (Lévy) |
| frequency of $a = 1, 2, 3$ | 0.413, 0.170, 0.093 | 0.415, 0.170, 0.093 (Gauss–Kuzmin) |

So, as far as one can see, $\ln C$ behaves like a generic real number. Its two features are:

- **the opening $[0; 1, 1, 1, 1, 1, 1, \dots]$**: about 0.5% of real numbers start this way, and it is the whole origin of the golden ratio flavour;
- **a giant quotient $a_{4091} = 73\,852$** (the next largest is $a_{136} = 8\,388$). By the theory above, $q_{4090}$, a number with about 2 100 digits, is an extraordinarily sharp near-period. It fails with density about $2/(73\,852\, q_{4090})$. Under Gauss–Kuzmin, a quotient this large appears among the first 11 000 terms with probability about 20%: rare, but not suspicious.
