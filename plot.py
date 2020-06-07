import matplotlib.pyplot as plt
from sobol import Sobol
import sys

if __name__ == "__main__":

    sobol = Sobol.load(sys.argv[1])
    samples = int(sys.argv[2])
    output = sys.argv[3]

    sobol_x = sobol.matrix(1, 32, 32, reverse=True)
    sobol_y = sobol.matrix(2, 32, 32, reverse=True)

    xs = map(lambda x: x / (2 ** 32), Sobol.generate(sobol_x))
    ys = map(lambda y: y / (2 ** 32), Sobol.generate(sobol_y))

    x = [next(xs) for _ in range(samples)]
    y = [next(ys) for _ in range(samples)]

    plt.scatter(x, y)
    plt.savefig(output)
