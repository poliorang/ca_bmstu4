from numpy import arange
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from numpy import linspace


class Point:
    def __init__(self, x=0, y=0, z=0, weight=1):
        self.x = x
        self.y = y
        self.z = z
        self.weight = weight

    def __str__(self):
        return f"|{self.x:^10.2f} | {self.y:^10.2f} | {self.z:^10.2f} | {self.weight:^10.2f} |"


def print_table(table):
    print(" __________________________________________________")
    print("|      X    |      Y     |     Z      |   weight   |")
    print("|--------------------------------------------------|")
    for i in range(len(table)):
        print(table[i])
    print("|__________________________________________________|")


def read_from_file(file_input):
    dots = list()
    with open(file_input, "r") as f:
        line = f.readline()
        while line:
            x, y, z, weight = map(float, line.split())
            dots.append(Point(x, y, z, weight))
            line = f.readline()
    return dots


def find_slae_matrix_3d(dots, n):
    rn = 0
    for i in range(0, n + 1):
        rn += i + 1
    res = [[0 for i in range(0, rn)] for j in range(0, rn)]
    col = [0 for i in range(0, rn)]
    v = 0
    for i in range(0, n + 1):
        for j in range(0, i + 1):
            for h in range(len(dots)):
                c = 0
                for k in range(0, n + 1):
                    for l in range(0, k + 1):
                        coef = dots[h].weight * dots[h].x ** (i - j) * dots[h].y ** j
                        # print(coef, sep=' ', end='')
                        res[v][c] += coef * dots[h].x ** (k - l) * dots[h].y ** l
                        c += 1
                # print('\n')
                col[v] += coef * dots[h].z
            v += 1

    for i in range(len(col)):
        res[i].append(col[i])
    return res


# Функция метод Гаусс
def method_gauss(matrix):
    n = len(matrix)
    for k in range(n):
        for i in range(k + 1, n):
            coeff = -(matrix[i][k] / matrix[k][k])
            for j in range(k, n + 1):
                matrix[i][j] += coeff * matrix[k][j]
    a = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            matrix[i][n] -= a[j] * matrix[i][j]
        a[i] = matrix[i][n] / matrix[i][i]
    return a


def get_polynomial_coefficients(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i == j:
                continue
            multiplication = matrix[j][i] / matrix[i][i]
            for k in range(0, len(matrix) + 1):
                matrix[j][k] -= multiplication * matrix[i][k]

    for i in range(len(matrix)):
        multiplication = matrix[i][i]
        for j in range(len(matrix[i])):
            matrix[i][j] /= multiplication

    return [matrix[i][-1] for i in range(len(matrix))]


def add_plot(coeffs, degree, label, subpl):
    my_x = []
    my_y = []
    my_z = []

    print(coeffs)
    for elx in x_arr:
        for ely in y_arr:
            elz = 0
            i = 0
            for k in range(0, degree + 1):
                for l in range(0, k + 1):
                    elz += coeffs[i] * elx ** (k - l) * ely ** l
                    i += 1
            my_x.append(elx)
            my_y.append(ely)
            my_z.append(elz)

    print(my_x)
    plt.plot(my_x, my_y, my_z, label=label)


def create_figure(table):
    global x_arr, y_arr
    fig = plt.figure()
    subpl = fig.add_subplot(111, projection='3d')

    table_x = [table[i].x for i in range(len(table))]
    table_y = [table[i].y for i in range(len(table))]
    table_z = [table[i].z for i in range(len(table))]

    subpl.scatter(table_x, table_y, table_z)

    x_arr = list(linspace(min(table_x), max(table_x), 5))
    y_arr = list(linspace(min(table_y), max(table_y), 5))

    return subpl


def add_table(table, label):
    table_x = [table[i].x for i in range(len(table))]
    table_y = [table[i].y for i in range(len(table))]
    table_z = [table[i].z for i in range(len(table))]

    # subpl.plot(table_x, table_y, table_z, 'o', label=label)


def draw_result():
    plt.legend()

    plt.xlabel('X') 
    plt.ylabel('Y')
    plt.ylabel('Z')

    plt.grid()
    plt.show()


if __name__ == "__main__":
    # filenames = input("Enter filenames: ").split()
    # label = input("Enter labels: ").split(',')

    filename = 'input.txt'
    degrees = int(input("Enter polynomial degree: "))

    points = read_from_file(filename)
    add_table(points, "Table")

    print_table(points)

    x_arr = []
    y_arr = []

    slae_matrix = find_slae_matrix_3d(points, degrees)
    coeffs = method_gauss(slae_matrix)
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    table_x = [points[i].x for i in range(len(points))]
    table_y = [points[i].y for i in range(len(points))]
    table_z = [points[i].z for i in range(len(points))]

    ax.scatter(table_x, table_y, table_z)

    x = list(linspace(min(table_x), max(table_x), 20))
    y = list(linspace(min(table_y), max(table_y), 20))
    res_x = []
    res_y = []
    res_z = []
    for elx in x:
        for ely in y:
            elz = 0
            i = 0
            for k in range(0, degrees + 1):
                for l in range(0, k + 1):
                    elz += coeffs[i] * elx ** (k - l) * ely ** l
                    i += 1
            res_x.append(elx)
            res_y.append(ely)
            res_z.append(elz)
    ax.plot_trisurf(res_x, res_y, res_z, color='pink', alpha=0.4)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    show()
