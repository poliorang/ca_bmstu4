
def print_data(data):
    print('   x   |   y  \n'
          '----------------')
    for point in data:
        print("%.3f  |  %.3f" % (point[0], point[1]))
    print('----------------\n')


def read_from_file(filename):
    try:
        with open(filename, "r") as f:
            data = [list(map(float, string[:-1].split())) for string in f.readlines()]
        print_data(data)
        return data
    except:
        return []


def newton(points, x, n):
    args = [p[0] for p in points]
    vals = [p[1] for p in points]

    difs = [vals[0]]
    col_num = len(vals)
    for i in range(1, col_num):
        for j in range(col_num - 1):
            vals[j] = (vals[j] - vals[j + 1]) / (args[j] - args[j + i])
        difs.append(vals[0])
        col_num -= 1

    result = 0
    multiplier = 1
    for i in range(n):
        result += (difs[i] * multiplier)
        multiplier *= (x - args[i])
    return result


def input_x():
    try:
        x = float(input("Введите х: "))
        return x
    except ValueError:
        print('Неверный ввод')
        exit(1)


def start_ratios(x, data, c, ksi, eta):
    print('Граничные условия:\n'
          '1: x0 = 0\n'         # вторая производная сплайна в левом краевом узле равна 0\n'
          '   xn = 0\n'         # вторая производная сплайна в правом краевом узле равна 0\n'
          '2: x0 = P3``(x0)\n'  # вторая производная сплайна в левом краевом узле равна\n'
                                # второй производной полинома Ньютона третьей степени\n'
          '   xn = 0\n'         # вторая производная сплайна в правом краевом узле равна 0\n'
          '3: x0 = P3``(x0)\n'  # вторая производная сплайна в левом краевом узле равна\n'
                                # второй производной полинома Ньютона третьей степени\n'
          '   xn = P3''(xn)')   # вторая производная сплайна в правом краевом узле равна\n'
                                # второй производной полинома Ньютона третьей степени\n'
    try:
        choice = int(input('Введите граничные условия: '))

        if choice == 1:
            return c, ksi, eta

        elif choice == 2:
            c[1] = newton(data, x, 3)
            ksi[2] = 1
            return c, ksi, eta

        elif choice == 3:
            c[-1] = newton(data, x, 3)
            c[1] = newton(data, x, 3)
            ksi[2] = 1
            return c, ksi, eta

        else:
            print('Неверный ввод')
            exit(1)

    except ValueError:
        print('Неверный ввод')
        exit(1)



# прямой ход нужен, чтобы вычислить
# прогоночные коэффициенты кси и эта
# для метода гаусса
# по дефолту кси2 и эта2 = 0 (из условия с1 = 0,
# в узлах значения многочлена и интерполируемой функции совпадают)
def straight_walk(ksi, eta, length, data, h):
    for i in range(3, length):
        f = -3 * ((data[i - 1][1] - data[i - 2][1]) / h[i - 1]
                  - (data[i - 2][1] - data[i - 3][1]) / h[i - 2])
        denominator = -2 * (h[i - 1] + h[i - 2]) - h[i - 2] * ksi[i - 1]
        ksi[i] = h[i - 1] / denominator
        eta[i] = (f + h[i - 2] * eta[i - 1]) / denominator


# обратный ход нужен, чтобы найти коэффициенты с
# они зависят от эты и кси6 которые мы нашли в прямом ходе
def forward_walk(ksi, eta, c):
    c[-2] = eta[-1]
    for i in range(len(c) - 2, 1, -1):
        c[i] = ksi[i + 1] * c[i + 1] + eta[i + 1]


# поиск коэффициентов d и b
def find_ratios(d, b, a, c, h):
    for i in range(1, len(d) - 1):
        d[i] = (c[i + 1] - c[i]) / 3 / h[i]
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] / 3 * (c[i + 1] + 2 * c[i])


# поиск интервала, в котором расположен введенный x
def find_section(points: list, x):
    if points[0][0] - points[1][0] < 0:
        i = 0
        while x > points[i][0] and i < len(points) - 1:
            i += 1
        return i
    else:
        i = 0
        while x < points[i][0] and i < len(points) - 1:
            i += 1
        return i


# поиск результата по формуле кубического полинома (xl - нижняя граница интервала)
def find_result(x, a, b, c, d, xl):
    return a + b * (x - xl) + c * (x - xl) ** 2 + d * (x - xl) ** 3
