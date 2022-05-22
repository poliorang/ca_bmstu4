from matplotlib import pyplot as plt

GAMMA = 0.2 * 10 ** (-2)
EPS = 10 ** (-4)
MAX_ITER = 10000
DELTA = 1.5 * 10 ** 3
A1 = 0.0134
B1 = 1
C1 = 4.35 * 10 ** (-4)
M1 = 1
R = 0.5
T0 = 300


def init_functions(k, a, p, f, i, T, a0, R):
    k[i] = A1 * (B1 + C1 * (T[i] ** M1))
    a[i] = a0 * (T[i] / DELTA - 1) ** 4 + GAMMA
    p[i] = 2 / R * a[i]
    f[i] = 2 * T0 / R * a[i]


# находим коэффициенты при левом краевом условии
def left_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, F0, a0, h, i=0):
    k[i + 1] = A1 * (B1 + C1 * (T[i + 1] ** M1))
    A[i] = 0
    B[i] = (k[i] + k[i + 1]) / 2 + p[i] * h ** 2
    C[i] = (k[i] + k[i + 1]) / 2
    D[i] = f[i] * h ** 2 + F0 * h
    A_[i] = 0
    B_[i] = B[i] + (A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1)) + 8 * a0 / (DELTA * R) * (
            T[i] / DELTA - 1) ** 3 * h ** 2) * T[i] - (A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1))) * T[
                i + 1] - (
                    8 * T0 * a0 / (DELTA * R) * (T[i] / DELTA - 1) ** 3 * h ** 2)
    C_[i] = -(A1 * C1 * M1 / 2 * (T[i + 1] ** (M1 - 1))) * T[i] + (A1 * C1 * M1 / 2 * (
            T[i + 1] ** (M1 - 1))) * T[i + 1] + C[i]
    D_[i] = - B[i] * T[i] + C[i] * T[i + 1] + D[i]


# находим коэффициенты при правом краевом условии
def right_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, a, a0, h, i):
    A[i] = (k[i - 1] + k[i]) / 2
    B[i] = (k[i - 1] + k[i]) / 2 + p[i] * h ** 2 + a[i] * h
    C[i] = 0
    D[i] = f[i] * h ** 2 + a[i] * T0 * h
    A_[i] = (A1 * C1 / 2 * M1 * (T[i - 1] ** (M1 - 1))) * T[i - 1] + A[i] - (A1 * C1 / 2 * M1 * (
            T[i - 1] ** (M1 - 1))) * T[i]
    B_[i] = -A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1)) * T[i - 1] + B[i] + (
            A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1)) + 8 * a0 / (DELTA * R) * (
            T[i] / DELTA - 1) ** 3 * h ** 2 + 4 * a0 / DELTA * (T[i] / DELTA - 1) ** 3 * h) * T[i] - (
                    8 * T0 * a0 / (DELTA * R) * (T[i] / DELTA - 1) ** 3 * h ** 2 + 4 * a0 / DELTA * (
                    T[i] / DELTA - 1) ** 3 * T0 * h)
    C_[i] = 0
    D_[i] = A[i] * T[i - 1] - B[i] * T[i] + D[i]


# находим коэффициенты
def non_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, a0, h, i):
    k[i + 1] = A1 * (B1 + C1 * (T[i + 1] ** M1))
    A[i] = (k[i - 1] + k[i]) / 2
    B[i] = (k[i - 1] + k[i]) / 2 + (k[i] + k[i + 1]) / 2 + p[i] * h ** 2
    C[i] = (k[i] + k[i + 1]) / 2
    D[i] = f[i] * h ** 2
    A_[i] = (A1 * C1 / 2 * M1 * (T[i - 1] ** (M1 - 1))) * T[i - 1] + A[i] - (
            A1 * C1 / 2 * M1 * (T[i - 1] ** (M1 - 1))) * T[i]
    B_[i] = -A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1)) * T[i - 1] + B[i] + (
            A1 * C1 * M1 * (T[i] ** (M1 - 1)) + 8 * a0 / (DELTA * R) * (
            T[i] / DELTA - 1) ** 3 * h ** 2) * T[i] - (
                    A1 * C1 * M1 / 2 * (T[i] ** (M1 - 1)) * T[i + 1]) - (
                    8 * T0 * a0 / (DELTA * R) * (T[i] / DELTA - 1) ** 3 * h ** 2)
    C_[i] = -(A1 * C1 * M1 / 2 * (T[i + 1] ** (M1 - 1))) * T[i] + (
            A1 * C1 * M1 / 2 * (T[i + 1] ** (M1 - 1))) * T[i + 1] + C[i]
    D_[i] = A[i] * T[i - 1] - B[i] * T[i] + C[i] * T[i + 1] + D[i]


# прогоночные коэффициенты, изначальные кси и эта равны нулю
def run_through_coeffs(ksi, eta, A_, B_, C_, D_, i):
    ksi[i] = C_[i - 1] / (B_[i - 1] - A_[i - 1] * ksi[i - 1])
    eta[i] = (A_[i - 1] * eta[i - 1] + D_[i - 1]) / (B_[i - 1] - A_[i - 1] * ksi[i - 1])


# обратный ход - вычисляем дельта y - расстояние между точками в сетке
def reverse_course(ksi, eta, N, dy, T, A_, B_, D_):
    max_dy = 0
    for i in range(N, -1, -1):
        if i < N:
            dy[i] = ksi[i + 1] * dy[i + 1] + eta[i + 1]
        else:
            dy[i] = (A_[i] * eta[i] + D_[i]) / (B_[i] - A_[i] * ksi[i])
        max_dy = max(max_dy, abs(dy[i] / T[i]))

    return max_dy


def show(N, h, a0, F0, T):
    x = [0 for j in range(N + 1)]
    for i in range(1, N + 1):
        x[i] = x[i - 1] + h
    fig, ax = plt.subplots()
    text = "N=" + str(N) + ", F0=" + str(F0) + ", a0=" + str(a0) + ", EPS=" + str(EPS)
    plt.title(text)
    ax.plot(x, T, '-')
    ax.set_xlabel('Расстояние от левого торца стержня, см')
    ax.set_ylabel('Температура, К')
    ax.grid(True)
    plt.show()


def main():
    while True:
        l = 10
        F0 = 50
        a0 = 1.94 * 10 ** (-2)
        N = 100
        print("1. F0 = 50 (подача тепла)\n"
              "2. F0 = -10 (отвод тепла)\n"
              "3. a0 = 3*a0 (ускорение подачи тепла)\n"
              "4. F0 = 0 (температура стержня не меняется)\n"
              "5. Собственные параметры\n"
              "0. Выход")
        test = int(input(">> "))
        if test == 0:
            print("Выход.")
            return
        elif test == 1:
            F0 = 50
        elif test == 2:
            F0 = -10
        elif test == 3:
            a0 = 3 * a0
        elif test == 4:
            F0 = 0
        elif test == 5:
            F0 = float(input("Введите F0: "))
            N = int(input("Введите число точек: "))
            coeff = float(input("Коэффициент, на который умножить a0 (a0 = 0.0194): "))
            a0 = a0 * coeff
        else:
            print("Данного пункта не существует.")
            return

        h = l / N
        T = [T0 for i in range(N + 1)]
        k = [0 for i in range(N + 1)]
        p = [0 for i in range(N + 1)]
        A = [0 for i in range(N + 1)]
        B = [0 for i in range(N + 1)]
        C = [0 for i in range(N + 1)]
        D = [0 for i in range(N + 1)]
        a = [0 for i in range(N + 1)]
        f = [0 for i in range(N + 1)]
        A_ = [0 for i in range(N + 1)]
        B_ = [0 for i in range(N + 1)]
        C_ = [0 for i in range(N + 1)]
        D_ = [0 for i in range(N + 1)]
        ksi = [0 for i in range(N + 1)]
        eta = [0 for i in range(N + 1)]
        dy = [0 for i in range(N + 1)]

        j = 0
        while j < MAX_ITER:
            for i in range(N + 1):

                init_functions(k, a, p, f, i, T, a0, R)  # задаем значения функций

                # получаем коэффициенты СЛАУ с трехдиагональной матрицей
                if i == 0:
                    left_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, F0, a0, h)
                elif i == N:
                    right_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, a, a0, h, N)
                else:
                    non_border_coeffs(k, A, B, C, D, A_, B_, C_, D_, T, f, p, a0, h, i)

                # прямой ход - находим кси и ета для метода Гаусса
                if i > 0:
                    run_through_coeffs(ksi, eta, A_, B_, C_, D_, i)

            # максимальная разность между значением и узлом
            max_dy = reverse_course(ksi, eta, N, dy, T, A_, B_, D_)

            if max_dy < EPS:
                break
            else:
                for i in range(N + 1):
                    T[i] += dy[i]

            j += 1

        show(N, h, a0, F0, T)


if __name__ == "__main__":
    main()
