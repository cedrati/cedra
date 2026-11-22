"""
Biquadratic equation solver: a x^4 + b x^2 + c = 0.

Author: Olivier Cedrati
License: CC BY-NC 4.0
Repository: https://github.com/cedrati/cedra

Side utility of the repository (C itself is a root of a quartic, but this
solver is general). Coefficients are normalised first, so tolerances are
relative and the result does not depend on the scale of (a, b, c).
"""

import math

REL_EPS = 1e-14


class CedraSolver:
    """Real roots of a x^4 + b x^2 + c = 0 through y = x^2."""

    def __init__(self):
        self.solutions_count = 0
        self.special_cases_count = 0

    def _roots_in_y(self, a, b, c):
        """Real roots of a y^2 + b y + c = 0 (coefficients already normalised)."""
        if a == 0:
            if b == 0:
                return []                      # c != 0 after normalisation
            self.special_cases_count += 1
            return [-c / b]

        discriminant = b * b - 4 * a * c
        if abs(discriminant) <= REL_EPS * (b * b + abs(4 * a * c)):
            self.special_cases_count += 1      # double root
            return [-b / (2 * a)]
        if discriminant < 0:
            return []

        # Numerically stable form: q = -(b + sign(b) sqrt(D)) / 2, roots q/a and c/q
        q = -(b + math.copysign(math.sqrt(discriminant), b)) / 2
        return [q / a, c / q]

    def solve(self, a, b, c):
        """
        Solve a x^4 + b x^2 + c = 0.

        Returns the distinct real solutions sorted in ascending order.
        Raises ValueError when a = b = c = 0 (every x is a solution).
        """
        self.solutions_count += 1

        if not all(math.isfinite(v) for v in (a, b, c)):
            return []
        scale = max(abs(a), abs(b), abs(c))
        if scale == 0:
            raise ValueError("0 = 0: every real number is a solution")
        a, b, c = a / scale, b / scale, c / scale

        solutions = []
        for y in self._roots_in_y(a, b, c):
            if abs(y) <= REL_EPS:
                solutions.append(0.0)
            elif y > 0:
                root = math.sqrt(y)
                solutions.extend([-root, root])

        unique = []
        for x in sorted(solutions):
            if not unique or abs(x - unique[-1]) > REL_EPS * max(1.0, abs(x)):
                unique.append(x)
        return unique

    def get_stats(self):
        """Return solver statistics."""
        return {
            'total_solutions': self.solutions_count,
            'special_cases': self.special_cases_count,
            'special_case_rate': self.special_cases_count / max(1, self.solutions_count) * 100,
        }


if __name__ == "__main__":
    solver = CedraSolver()

    test_cases = [
        (1, -5, 4),            # +-1, +-2
        (1, -2, 1),            # +-1 (double root in y)
        (1, 0, -1),            # +-1
        (1, -10, 9),           # +-1, +-3
        (2, -8, 6),            # +-1, +-sqrt(3)
        (1, 1, 1),             # no real solution
        (1, -3, 0),            # 0, +-sqrt(3)
        (1e-20, -5e-20, 4e-20),  # same as the first one, tiny scale
    ]

    print("Biquadratic solver test results")
    print("=" * 50)
    for i, (a, b, c) in enumerate(test_cases, 1):
        solutions = solver.solve(a, b, c)
        print(f"Test {i}: {a}x^4 + {b}x^2 + {c} = 0")
        if solutions:
            print(f"  Solutions: {[f'{x:.6f}' for x in solutions]}")
            scale = max(abs(a), abs(b), abs(c))
            max_error = max(abs(a * x**4 + b * x**2 + c) / scale for x in solutions)
            print(f"  Max relative residual: {max_error:.2e}")
        else:
            print("  No real solutions")

    stats = solver.get_stats()
    print("\nSolver statistics:")
    print(f"  Total equations solved: {stats['total_solutions']}")
    print(f"  Special cases detected: {stats['special_cases']}")
    print(f"  Special case rate: {stats['special_case_rate']:.1f}%")
