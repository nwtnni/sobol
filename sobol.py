import copy
from functools import reduce
import sys


class Sobol:

    # Compute the Sobol generator matrix V:
    #
    #                          A
    #          V            +-----+
    # +-----+-----+-----+   | a_1 |
    # | v_1 | ... | v_n | * | ... |
    # +-----+-----+-----+   | a_n |
    #                       +-----+
    #
    # Here `v_1, ..., v_n` are `m`-bit numbers.
    #
    def matrix(self, dimension, m, n, reverse=True, binary=True):
        v = self.invert(self.directions(dimension, n), m)

        if binary and reverse:
            return map(lambda v_i: "{:0{}b}".format(v_i, m)[::-1], v)
        elif binary:
            return map(lambda v_i: "{:0{}b}".format(v_i, m), v)
        elif reverse:
            v = map(lambda v_i: int("{:0{}b}".format(v_i, m)[::-1], 2), v)
            return map(lambda v_i: "{:0{}x}".format(v_i, (m + 3) // 4), v)
        else:
            return map(lambda v_i: "{:0{}x}".format(v_i, (m + 3) // 4), v)

    # Compute inverse direction numbers:
    #
    # v_1 = m_1 / 2^1
    # v_2 = m_2 / 2^2
    # ...
    # v_n = m_n / 2^n
    #
    def invert(self, m, precision):
        return [
            int("{:0{}b}".format(m_i, i)[:precision][::-1], 2)
            for (i, m_i) in enumerate(m, 1)
        ]

    # Compute the first `n` direction numbers `m_1, ..., m_n` for `dimension`
    def directions(self, dimension, n):
        if dimension == 0:
            return [1 for _ in range(n)]

        s = self.s[dimension - 1]
        a = self.a[dimension - 1]
        m = copy.deepcopy(self.m_i[dimension - 1])

        # Compute `a XOR b`
        def xor(a, b):
            return a ^ b

        # Compute the 1-indexed `n`th bit of integer `a`
        def bit(a, n):
            return (a >> (n - 1)) & 1

        while len(m) < n:
            # m_i = 2^1 * a_1 * m_{i - 1}
            #     ^ 2^2 * a_2 * m_{i - 2}
            #     ^ ...
            #     ^ 2^{s - 1} * a_{s - 1} * m_{i - d + 1}
            m_i = reduce(xor, [bit(a, d) * m[-d] << d for d in range(1, s)], 0)

            #     ^ 2^s * m_{i - d}
            m_i ^= m[-s] << s

            #     ^ m_{i - d}
            m_i ^= m[-s]

            m.append(m_i)

        return m

    # Load a series of primitive polynomials and direction numbers.
    #
    # Requires input in the following format:
    #
    # ```txt
    # d    s    a    m_i
    # d_0  s_0  a_0  m_i_0_0 m_i_0_1 m_i_0_2 ...
    # d_1  s_1  a_1  m_i_1_0 m_i_1_1 m_i_1_2 ...
    # .    .    .    .
    # .    .    .    .
    # .    .    .    .
    # d_n  s_n  a_n  m_i_n_0 m_i_n_1 m_i_n_2 ...
    # ```
    @staticmethod
    def load(path):
        s = []
        a = []
        m_i = []

        with open(path, "r") as file:

            # Skip header
            next(file)

            for line in file:
                _, _s, _a, *_m_i = [int(i) for i in line.strip().split()]
                s.append(_s)
                a.append(_a)
                m_i.append(_m_i)

        return Sobol(s, a, m_i)

    def __init__(self, s, a, m_i):
        self.s = s
        self.a = a
        self.m_i = m_i


if __name__ == "__main__":
    # Test case from Bratley, Fox
    #
    # Polynomial: x^3 + x^2 + 1
    # Initial values: m_1 = 1, m_2 = 3, m_3 = 7
    #
    # Using notation from Joe, Kuo, we have:
    #
    # s = 3
    # a = 0b010
    # m_i = 1, 3, 7
    sobol = Sobol([3], [2], [[1, 3, 7]])
    assert(sobol.directions(1, 6) == [1, 3, 7, 5, 7, 43])

    sobol = Sobol.load(sys.argv[1])
    for i in range(4):
        for line in sobol.matrix(i, 32, 52, reverse=True, binary=False):
            print(line)
        print()
