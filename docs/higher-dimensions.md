# Cedra in higher dimensions: the Galois orbit

*Olivier Cedrati, November 2026*

A one-dimensional Kronecker sequence cannot fill the square: the pairs $(\{n\alpha\}, \{(n+1)\alpha\})$ lie on two lines. To sample in $d$ dimensions one needs $d$ numbers $\theta_1, \dots, \theta_d$ such that $1, \theta_1, \dots, \theta_d$ are linearly independent over $\mathbb{Q}$. Cedra comes with natural partners, its **Galois conjugates**. This note proves that their logarithms are independent, using a prime number hidden in the minimal polynomial, and then measures honestly how good the resulting point sets are. Numbers are reproduced by [`cedra_multidim.py`](../cedra_multidim.py) and covered by [`test_cedra.py`](../test_cedra.py).

## The Galois orbit

$C$ lives in $K = \mathbb{Q}(\sqrt{2}, \sqrt{3})$, whose Galois group has four elements (the sign changes of $\sqrt{2}$ and $\sqrt{3}$). The four conjugates are the roots of $P(x) = 4x^4 + 32x^3 + 36x^2 - 112x - 167$:

| | expression | value | $\ln\lvert C_i\rvert$ |
|---|---|---|---|
| $C_1 = C$ | $\sqrt{3} + \tfrac{3}{2}\sqrt{2} - 2$ | 1.853371 | 0.617006 |
| $C_2$ | $\sqrt{3} - \tfrac{3}{2}\sqrt{2} - 2$ | -2.389270 | 0.870988 |
| $C_3$ | $-\sqrt{3} + \tfrac{3}{2}\sqrt{2} - 2$ | -1.610730 | 0.476688 |
| $C_4 = -(C + 4)$ | $-\sqrt{3} - \tfrac{3}{2}\sqrt{2} - 2$ | -5.853371 | 1.767018 |

## Theorem D

$$\boxed{\ 1,\ \ln\lvert C_1\rvert,\ \ln\lvert C_2\rvert,\ \ln\lvert C_3\rvert,\ \ln\lvert C_4\rvert \text{ are linearly independent over } \mathbb{Q}\ }$$

Consequently, by Weyl's criterion, the sequence $\big(\{n \ln\lvert C_1\rvert\}, \dots, \{n \ln\lvert C_4\rvert\}\big)$ is equidistributed in $[0,1)^4$, and so is every sub-family in lower dimension.

### The prime 167

The key is the constant term of $P$.

- $2C$ is an algebraic integer, and its norm is $N(2C) = 2^4 \prod C_i = 16 \cdot (-167/4) = -668 = -2^2 \cdot 167$.
- $167$ is prime, $167 \equiv 7 \pmod 8$ and $167 \equiv -1 \pmod{12}$, so 2 and 3 are both squares modulo 167. Hence 167 **splits completely** in $K$, into four distinct prime ideals of norm 167, which the Galois group permutes simply transitively. Equivalently, $P$ has four distinct roots modulo 167: $0, 39, 124, 163$.

**Lemma.** Each $C_i$ is divisible by exactly one prime $\mathfrak{p}_i$ above 167, with multiplicity 1, and the $\mathfrak{p}_i$ are pairwise distinct.

*Proof.* The ideal $(2C_1)$ has norm $668 = 4 \cdot 167$, and each prime above 167 has norm 167. So exactly one of them, $\mathfrak{p}_1$, divides $2C_1$, with exponent 1; since 2 is prime to 167, the same holds for $C_1$. Applying the Galois automorphism $g_i$ that sends $C_1$ to $C_i$, the prime $\mathfrak{p}_i = g_i(\mathfrak{p}_1)$ is the only one above 167 that divides $C_i$. The $g_i$ are distinct and act simply transitively, so the $\mathfrak{p}_i$ are distinct. $\square$

### Proof of Theorem D

Suppose $c_0 + \sum_i c_i \ln\lvert C_i\rvert = 0$ with integers $c_i$ (clear denominators).

- If $c_0 \neq 0$: $e^{-c_0} = \prod_i \lvert C_i\rvert^{c_i}$ is algebraic, which contradicts the Lindemann theorem ($e^a$ is transcendental for every nonzero algebraic $a$).
- If $c_0 = 0$: $\prod_i C_i^{c_i} = \pm 1$. Taking the $\mathfrak{p}_j$-adic valuation, the lemma gives $c_j = 0$ for every $j$. $\square$

**A two-line proof in dimension 2.** For the pair $(C_1, C_2)$ alone, there is an elementary argument. The automorphism $\sqrt{2} \mapsto -\sqrt{2}$ swaps $C_1$ and $C_2$. So a relation $C_1^a C_2^b = \pm 1$ implies $C_2^a C_1^b = \pm 1$. Multiplying and dividing gives $(C_1C_2)^{a+b} = \pm 1$ and $(C_1/C_2)^{a-b} = \pm 1$. Since $\lvert C_1C_2\rvert = 4\sqrt{3} - 2.5 \approx 4.43$ and $\lvert C_1/C_2\rvert \approx 0.78$ are not 1, we get $a = b = 0$.

## How good are these point sets?

Theorem D guarantees convergence, not speed. We measure the $L^2$ star discrepancy (Warnock's exact formula; lower is better) against standard constructions.

**Dimension 2**

| point set | N = 500 | N = 2000 |
|---|---|---|
| Cedra $(\ln C_1, \ln\lvert C_2\rvert)$ | $7.3 \times 10^{-3}$ | $1.6 \times 10^{-3}$ |
| Cedra $(\ln C_1, \ln\lvert C_3\rvert)$ | $1.6 \times 10^{-2}$ | $1.9 \times 10^{-3}$ |
| R2 (plastic number) | $3.7 \times 10^{-3}$ | $\mathbf{6.6 \times 10^{-4}}$ |
| $(\sqrt{2}, \sqrt{3})$ | $3.9 \times 10^{-3}$ | $1.1 \times 10^{-3}$ |
| Halton (2, 3) | $\mathbf{2.6 \times 10^{-3}}$ | $7.4 \times 10^{-4}$ |
| pseudo-random | $1.8 \times 10^{-2}$ | $5.9 \times 10^{-3}$ |

**Dimension 4**, N = 1500

| point set | $L^2$ star discrepancy |
|---|---|
| Cedra, all four conjugates | $4.6 \times 10^{-3}$ |
| $(\sqrt{2}, \sqrt{3}, \sqrt{5}, \sqrt{7})$ | $2.8 \times 10^{-3}$ |
| Halton (2, 3, 5, 7) | $\mathbf{1.4 \times 10^{-3}}$ |
| pseudo-random | $5.2 \times 10^{-3}$ |

**Simultaneous approximation.** The 2D analogue of the near-periods is an $n$ for which both $\lVert n\theta_1\rVert$ and $\lVert n\theta_2\rVert$ are small. Dirichlet's theorem gives infinitely many $n$ with $n \max(\lVert n\theta_1\rVert, \lVert n\theta_2\rVert)^2 \le 1$. A pair is *badly approximable* when this quantity stays bounded below, which is what makes a 2D Kronecker sequence uniform. The smallest values up to $n = 2 \times 10^6$:

| pair | smallest $n \max(\cdot)^2$ | at $n$ |
|---|---|---|
| R2 (plastic number, cubic field) | 0.172 | 31 572 |
| Cedra $(\ln C_1, \ln\lvert C_2\rvert)$ | 0.076 | 248 |
| $(\sqrt{2}, \sqrt{3})$ | 0.015 | 326 491 |

R2 comes from a cubic field, so the pair is provably badly approximable (Perron, Cassels). That explains why it wins. The Cedra pair is transcendental, so, like almost every pair (Khintchine), it is expected *not* to be badly approximable, though this is not proven.

## Conclusion

- **Proven:** the Galois orbit of $C$ yields a 4-dimensional Kronecker sequence that is equidistributed, through a short argument combining the prime $167$ and Lindemann's theorem.
- **Measured:** these point sets beat random sampling in 2D, but not standard low-discrepancy constructions (R2, Halton). In 4D, at N = 1500, they are at the level of random sampling.
- For practical quasi-Monte Carlo in several dimensions, use R2 or Halton/Sobol. The interest of the Cedra orbit is structural, not numerical.
