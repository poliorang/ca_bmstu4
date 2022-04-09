from numpy import arange
import matplotlib.pyplot as plt


class Point:
    def __init__(self, x=0, y=0, weight=1):
        self.x = x
        self.y = y
        self.weight = weight

    def __str__(self):
        return f"|{self.x:^10.2f} | {self.y:^10.2f} | {self.weight:^10.2f} |"


def print_table(table):
    print(" _____________________________________")
    print("|      X    |      Y     |   weight   |")
    print("|-------------------------------------|")
    for i in range(len(table)):
        print(table[i])
    print("|_____________________________________|")


def read_from_file(file_input):
    dots = list()
    with open(file_input, "r") as f:
        line = f.readline()
        while line:
            x, y, weight = map(float, line.split())
            dots.append(Point(x, y, weight))
            line = f.readline()
    return dots


def append_right_side(matrix, dots):
    for i in range(len(matrix)):
        res = 0
        for j in range(len(dots)):
            res += dots[j].weight * dots[j].y * (dots[j].x ** i)
        matrix[i].append(res)


def get_coefficient(dots, degree):
    coefficient = 0
    for i in range(len(dots)):
        coefficient += dots[i].weight * (dots[i].x ** degree)
    return coefficient


def find_slae_matrix_3d(dots, degree):
    matrix = []
    for h in range(len(dots)):
        coefficient = 0
        row = []
        for i in range(degree + 1):
            for j in range(degree + 1):
                coefficient += dots[i].weight * (dots[i].y ** (i - j)) * (dots[i].x ** degree)
            row.append(coefficient)
        matrix.append(row)
    return matrix


def find_slae_matrix(dots, degree):
    matrix = [[get_coefficient(dots, j + i)  # получится degree * 2
              for i in range(degree + 1)]
              for j in range(degree + 1)]
    append_right_side(matrix, dots)
    return matrix


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


def add_plot(coeffs, label, start, end):
    my_x = list()
    my_y = list()
    step = (end - start) / 1000
    for x in arange(start, end + step, step):
        my_x.append(x)
        y = 0
        for i in range(len(coeffs)):
            y += coeffs[i] * x ** i
        my_y.append(y)

    plt.plot(my_x, my_y, label=label)


def add_table(table, label):
    table_x = [table[i].x for i in range(len(table))]
    table_y = [table[i].y for i in range(len(table))]

    plt.plot(table_x, table_y, 'o', label=label)


def draw_result():
    plt.legend()

    plt.xlabel('X') 
    plt.ylabel('Y')

    plt.grid()
    plt.show()


if __name__ == "__main__":
    # filenames = input("Enter filenames: ").split()
    # label = input("Enter labels: ").split(',')

    filename = 'input2.txt'
    degrees = list(map(int, input("Enter polynomial degree: ").split()))

    points = read_from_file(filename)
    add_table(points, "Table")

    print_table(points)

    for j in range(len(degrees)):
        slae_matrix = find_slae_matrix(points, degrees[j])
        for i in range(len(slae_matrix[0])):
            print(*slae_matrix[i])
        coeffs = get_polynomial_coefficients(slae_matrix)
        print(*slae_matrix)
        print(coeffs)
        add_plot(coeffs, f"n = {degrees[j]}",
                 points[0].x, points[-1].x)

    draw_result()
