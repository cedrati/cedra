"""Unit tests: python -m unittest test_cedra"""

from decimal import Decimal
import random
import unittest

from cedra import (ALPHA_DEC, CEDRA_DEC, INV_PHI_DEC, MINIMAL_POLYNOMIAL,
                   continued_fraction, convergent_denominators, distance_to_integer,
                   eval_poly, factor_complexity, frac, kronecker, mismatch_density,
                   mismatch_positions, model_set, star_discrepancy, window_word)
from cedra_biquadratic import CedraSolver


class TestConstant(unittest.TestCase):
    def test_value(self):
        self.assertAlmostEqual(float(CEDRA_DEC), 1.8533711511285199, places=15)

    def test_three_terms(self):
        simplified = Decimal(3).sqrt() + Decimal(3) / 2 * Decimal(2).sqrt() - 2
        self.assertLess(abs(CEDRA_DEC - simplified), Decimal("1e-55"))

    def test_minimal_polynomial(self):
        self.assertLess(abs(eval_poly(MINIMAL_POLYNOMIAL, CEDRA_DEC)), Decimal("1e-50"))

    def test_delta_identity_is_tautological(self):
        # 2 sqrt(3) + 3 sqrt(2) - 4 == 2C, hence C * delta == phi for any C
        lhs = 2 * Decimal(3).sqrt() + 3 * Decimal(2).sqrt() - 4
        self.assertLess(abs(lhs - 2 * CEDRA_DEC), Decimal("1e-55"))

    def test_continued_fraction(self):
        self.assertEqual(continued_fraction(ALPHA_DEC, 9), [0, 1, 1, 1, 1, 1, 1, 3, 29])


class TestQuasicrystal(unittest.TestCase):
    N = 5000

    def test_frac_negative(self):
        self.assertEqual(frac(Decimal("-0.25")), Decimal("0.75"))

    def test_model_set_complexity_2n(self):
        word = window_word(self.N, INV_PHI_DEC)
        for k in (5, 20, 60):
            self.assertEqual(factor_complexity(word, k), 2 * k)

    def test_sturmian_complexity_n_plus_1(self):
        word = window_word(self.N, ALPHA_DEC)
        for k in (5, 20, 60):
            self.assertEqual(factor_complexity(word, k), k + 1)

    def test_convergent_denominators(self):
        self.assertEqual(convergent_denominators(ALPHA_DEC, 11),
                         [1, 1, 2, 3, 5, 8, 13, 47, 1376, 1423, 5645])

    def test_29_over_47_is_a_convergent(self):
        self.assertLess(abs(ALPHA_DEC - Decimal(29) / 47), Decimal(1) / (47 * 1376))

    def test_near_period_mismatch_rate(self):
        # Shifting by q changes a bit only near the two window edges: rate = 2 ||q alpha||
        word = window_word(50000, INV_PHI_DEC)
        for q in (47, 1423):
            rate = sum(word[i] != word[i + q] for i in range(len(word) - q)) / (len(word) - q)
            self.assertAlmostEqual(rate / (2 * float(distance_to_integer(q * ALPHA_DEC))), 1, delta=0.1)

    def test_gaps(self):
        pts = model_set(self.N)
        self.assertEqual({b - a for a, b in zip(pts, pts[1:])}, {1, 2})


class TestNearPeriods(unittest.TestCase):
    N = 100000

    def test_theorem_a(self):
        for window in (INV_PHI_DEC, ALPHA_DEC):
            word = window_word(self.N, window)
            for q in (3, 4, 13, 47, 100, 1376):
                measured = len(mismatch_positions(word, q)) / (self.N - q)
                self.assertAlmostEqual(measured / float(mismatch_density(q, window)), 1, delta=0.1)

    def test_theorem_b_sturmian_first_mismatch(self):
        word = window_word(10000, ALPHA_DEC)
        # (q_k, first mismatch): q_{k+1} when q_k alpha - p_k > 0, else 1
        for q, first in ((5, 8), (8, 1), (13, 47), (47, 1), (1376, 1423), (1423, 1)):
            self.assertEqual(mismatch_positions(word, q)[0], first)


class TestDiscrepancy(unittest.TestCase):
    def test_golden_ratio_beats_cedra(self):
        n = 10000
        self.assertLess(star_discrepancy(kronecker(n, INV_PHI_DEC)),
                        star_discrepancy(kronecker(n, ALPHA_DEC)))

    def test_known_value(self):
        self.assertAlmostEqual(star_discrepancy([0.5]), 0.5)


class TestSolver(unittest.TestCase):
    def test_scale_invariance(self):
        s = CedraSolver()
        for scale in (1e-20, 1.0, 1e20):
            got = s.solve(scale, -5 * scale, 4 * scale)
            self.assertEqual(len(got), 4)
            for x, expected in zip(got, (-2, -1, 1, 2)):
                self.assertAlmostEqual(x, expected, places=12)

    def test_special_case_counter(self):
        s = CedraSolver()
        s.solve(1, -5, 4)
        self.assertEqual(s.get_stats()['special_cases'], 0)
        s.solve(1, -2, 1)
        self.assertEqual(s.get_stats()['special_cases'], 1)

    def test_zero_root(self):
        self.assertEqual(CedraSolver().solve(1, -4, 0), [-2.0, 0.0, 2.0])

    def test_all_zero(self):
        with self.assertRaises(ValueError):
            CedraSolver().solve(0, 0, 0)

    def test_random_against_residual(self):
        s, rng = CedraSolver(), random.Random(0)
        for _ in range(2000):
            a, b, c = (rng.uniform(-10, 10) for _ in range(3))
            for x in s.solve(a, b, c):
                self.assertLess(abs(a * x**4 + b * x**2 + c), 1e-9 * max(1, x**4))


if __name__ == "__main__":
    unittest.main()
