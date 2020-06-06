import sys


class Sobol:
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
    def __init__(self, path):
        self.s = []
        self.a = []
        self.m_i = []

        with open(path, "r") as file:

            # Skip header
            next(file)

            for line in file:
                _, _s, _a, *_m_i = [int(i) for i in line.strip().split()]
                self.s.append(_s)
                self.a.append(_a)
                self.m_i.append(_m_i)


if __name__ == "__main__":
    path = sys.argv[1]
    sobol = Sobol(path)
