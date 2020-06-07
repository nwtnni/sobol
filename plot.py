import matplotlib.pyplot as plt
from sobol import Sobol
from random import random
import sys


def random_generator():
    while True:
        yield random()


def sobol_generator(sobol, dimension):
    matrix = sobol.matrix(dimension, 32, 32, reverse=True)
    return map(lambda x: x / (2 ** 32), Sobol.generate(matrix))


if __name__ == "__main__":

    samples = int(sys.argv[1])
    uniform = sys.argv[2] == "true"
    output = None

    if uniform:
        output = "uniform-"
    else:
        output = "sobol-"
    output += str(samples) + ".png"

    xg = None
    yg = None

    if uniform:
        xg = random_generator()
        yg = random_generator()
    else:
        sobol = Sobol.load("data/new-joe-kuo-6.21201")
        xg = sobol_generator(sobol, 0)
        yg = sobol_generator(sobol, 1)

    xs = [next(xg) for _ in range(samples)]
    ys = [next(yg) for _ in range(samples)]

    plt.scatter(xs, ys)
    plt.savefig(output)
