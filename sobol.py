import copy
from functools import reduce
import sys


class Sobol:

    # Returns a generator that yields sequential samples from `matrix`
    # starting at `index`.
    #
    # Faster than calling `.sample()` one at a time due to iterating
    # over samples in Gray code order.
    @staticmethod
    def generate(matrix, index=0):
        sample = Sobol.sample(matrix, gray(index))
        while True:
            index += 1
            sample ^= matrix[trailing_zeros(index)]
            yield sample

    # Returns sample `index` from `matrix`.
    #
    # Interpreting the output depends on how `matrix` was constructed.
    # Assuming output precision `o` and the flag `reverse = True`, then
    # the floating point value of the output is:
    #
    # ```python
    # float(output) / (2 ** o)
    # ```
    #
    # Otherwise, if `reverse = False`, the output directly represents
    # the radical inverse in base 2:
    #
    # ```txt
    # output_1 * 2^{-1} + output_2 * 2^{-2} + ... + output_o * 2^{-o}
    # ```
    #
    # Where `output_i` is the `i`th bit of `output` (counting left from
    # the least-significant bit). Here's another [explanation from PBR][pbr].
    #
    # [pbr]: http://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/(0,_2)-Sequence_Sampler.html#eq:generator-base-scale
    @staticmethod
    def sample(matrix, index):
        column = 0
        output = 0
        while index > 0:
            if index & 1 > 0:
                output ^= matrix[column]
            index >>= 1
            column += 1
        return output

    # Compute the Sobol generator matrix V for dimension `d`.
    #
    # ```txt
    #                          A
    #          V            +-----+
    # +-----+-----+-----+   | a_1 |
    # | v_1 | ... | v_n | * | ... |
    # +-----+-----+-----+   | a_n |
    #                       +-----+
    # ```
    #
    # The resulting matrix can produce up to `i` bits of output
    # for sample indices from `0` to `2^i - 1`. The output takes
    # on discrete values in steps of `1 / 2^o`.
    def matrix(self, d, i, o, reverse=False):
        m = self.directions(d, i)
        v = self.invert(m, o)
        if reverse:
            v = list(map(lambda v_i: reverse_bits(v_i, i), v))
        return v

    # Compute inverse direction numbers `v1, ... v_n` according to:
    #
    # ```txt
    # v_1 = m_1 / 2^1
    # v_2 = m_2 / 2^2
    # ...
    # v_n = m_n / 2^n
    # ```
    #
    # Only the `p` most significant bits of `m_1, ... m_n` will be used.
    def invert(self, m, p):
        return [reverse_bits(m_i, i, p) for (i, m_i) in enumerate(m, 1)]

    # Compute the first `n` direction numbers `m_1, ..., m_n` for dimension `d`
    def directions(self, d, n):
        if d == 0:
            return [1 for _ in range(n)]

        s = self.s[d - 1]
        a = self.a[d - 1]
        m = copy.deepcopy(self.m_i[d - 1])

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


# Compute the `n`th integer in [Gray code][gc] order.
#
# [gc]: https://en.wikipedia.org/wiki/Gray_code
def gray(n):
    return n ^ (n >> 1)


# Compute the number of trailing 0's in the binary representation of `n`.
def trailing_zeros(n):
    count = 0
    while n & 1 == 0:
        count += 1
        n >>= 1
    return count


# Take the `p` most significant bits of `n` and reverse them,
# adding leading zeros to pad to bit length `b` if necessary.
def reverse_bits(n, b, p=None):
    if p is None:
        p = b
    return int("{:0{}b}".format(n, b)[:p][::-1], 2)


if __name__ == "__main__":
    # Test cases from Bratley, Fox
    #
    # Polynomial: x^3 + x^2 + 1
    # Initial values: m_1 = 1, m_2 = 3, m_3 = 7
    #
    # Using notation from Joe, Kuo, we have:
    #
    # s = 3
    # a = 0b010
    # m_i = 1, 3, 7
    #
    sobol = Sobol([3], [2], [[1, 3, 7]])
    assert(sobol.directions(1, 6) == [1, 3, 7, 5, 7, 43])

    matrix = sobol.matrix(1, 5, 6)
    assert(Sobol.sample(matrix, gray(23)) == int("10001", 2))

    generator = Sobol.generate(matrix, 0)
    for (index, sample) in enumerate(generator, 1):
        if index == 1:
            assert(sample == int("01", 2))
        elif index == 2:
            assert(sample == int("10", 2))
        elif index == 3:
            assert(sample == int("11", 2))
        elif index == 23:
            assert(sample == int("10001", 2))
        elif index > 23:
            break

    sobol = Sobol.load(sys.argv[1])
    for i in range(4):
        for line in sobol.matrix(i, 32, 52, reverse=True):
            print("{:0{}b}".format(line, 32))
        print()
