from matplotlib import pyplot as plt


def main():
    EPS = 10 ** (-4)
    l = 10
    T0 = 300
    R = 0.5
    F0 = 50
    a1 = 0.0134
    b1 = 1
    c1 = 4.35 * 10 ** (-4)
    m1 = 1
    # a0 = -1.94 * 10**(-2)
    a0 = 1.94 * 10 ** (-2)
    delta = 1.5 * 10 ** 3
    gamma = 0.2 * 10 ** (-2)
    N = 100
    while True:
        F0 = 50
        a0 = 1.94 * 10 ** (-2)
        N = 100
        print("1) Тест 1 - F0 = 50 (подача тепла)\n"
              "2) Тест 2 - F0 = -10 (отвод тепла)\n"
              "3) Тест 3 - a0 = 3*a0 (ускорение подачи тепла)\n"
              "4) Тест 4 - F0 = 0 (температура стержня не меняется)\n"
              "5) Пользовательский режим\n"
              "0) Выход")
        test = int(input("Выберите тест (0-5): "))
        if test == 2:
            F0 = -10
        elif test == 3:
            a0 = 3 * a0
        elif test == 4:
            F0 = 0
        elif test == 5:
            F0 = float(input("Введите F0 (по умолч. F0 = 50): "))
            N = int(input("Введите число точек (по умолч. N = 100): "))
            coeff = float(input("Введите во сколько раз изменить a0 (по умолч. a0 = 0.0194): "))
            a0 = a0 * coeff
        elif test == 0:
            print("Выполнен выход!")
            return
        elif test == 1:
            F0 = 50
        else:
            print("Нет такого теста!")
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
        max_dy = 0
        max_iter = 10000
        j = 0

        while j < max_iter:
            for i in range(N + 1):
                k[i] = a1 * (b1 + c1 * (T[i] ** m1))
                a[i] = a0 * (T[i] / delta - 1) ** 4 + gamma
                p[i] = 2 / R * a[i]
                f[i] = 2 * T0 / R * a[i]
                if i == 0:
                    k[i + 1] = a1 * (b1 + c1 * (T[i + 1] ** m1))
                    A[i] = 0
                    B[i] = (k[i] + k[i + 1]) / 2 + p[i] * h ** 2
                    C[i] = (k[i] + k[i + 1]) / 2
                    D[i] = f[i] * h ** 2 + F0 * h
                    A_[i] = 0
                    B_[i] = B[i] + (a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1)) + 8 * a0 / (delta * R) * (
                            T[i] / delta - 1) ** 3 * h ** 2) * T[i] - (a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1))) * T[
                                i + 1] - (
                                    8 * T0 * a0 / (delta * R) * (T[i] / delta - 1) ** 3 * h ** 2)
                    C_[i] = -(a1 * c1 * m1 / 2 * (T[i + 1] ** (m1 - 1))) * T[i] + (a1 * c1 * m1 / 2 * (
                            T[i + 1] ** (m1 - 1))) * T[i + 1] + C[i]
                    D_[i] = - B[i] * T[i] + C[i] * T[i + 1] + D[i]
                elif i == N:
                    A[i] = (k[i - 1] + k[i]) / 2
                    B[i] = (k[i - 1] + k[i]) / 2 + p[i] * h ** 2 + a[i] * h
                    C[i] = 0
                    D[i] = f[i] * h ** 2 + a[i] * T0 * h
                    A_[i] = (a1 * c1 / 2 * m1 * (T[i - 1] ** (m1 - 1))) * T[i - 1] + A[i] - (a1 * c1 / 2 * m1 * (
                            T[i - 1] ** (m1 - 1))) * T[i]
                    B_[i] = -a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1)) * T[i - 1] + B[i] + (
                            a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1)) + 8 * a0 / (delta * R) * (
                            T[i] / delta - 1) ** 3 * h ** 2 + 4 * a0 / delta * (T[i] / delta - 1) ** 3 * h) * T[i] - (
                                    8 * T0 * a0 / (delta * R) * (T[i] / delta - 1) ** 3 * h ** 2 + 4 * a0 / delta * (
                                    T[i] / delta - 1) ** 3 * T0 * h)
                    C_[i] = 0
                    D_[i] = A[i] * T[i - 1] - B[i] * T[i] + D[i]
                else:
                    k[i + 1] = a1 * (b1 + c1 * (T[i + 1] ** m1))
                    A[i] = (k[i - 1] + k[i]) / 2
                    B[i] = (k[i - 1] + k[i]) / 2 + (k[i] + k[i + 1]) / 2 + p[i] * h ** 2
                    C[i] = (k[i] + k[i + 1]) / 2
                    D[i] = f[i] * h ** 2
                    A_[i] = (a1 * c1 / 2 * m1 * (T[i - 1] ** (m1 - 1))) * T[i - 1] + A[i] - (
                            a1 * c1 / 2 * m1 * (T[i - 1] ** (m1 - 1))) * T[i]
                    B_[i] = -a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1)) * T[i - 1] + B[i] + (
                            a1 * c1 * m1 * (T[i] ** (m1 - 1)) + 8 * a0 / (delta * R) * (
                            T[i] / delta - 1) ** 3 * h ** 2) * T[i] - (
                                        a1 * c1 * m1 / 2 * (T[i] ** (m1 - 1)) * T[i + 1]) - (
                                    8 * T0 * a0 / (delta * R) * (T[i] / delta - 1) ** 3 * h ** 2)
                    C_[i] = -(a1 * c1 * m1 / 2 * (T[i + 1] ** (m1 - 1))) * T[i] + (
                            a1 * c1 * m1 / 2 * (T[i + 1] ** (m1 - 1))) * T[i + 1] + C[i]
                    D_[i] = A[i] * T[i - 1] - B[i] * T[i] + C[i] * T[i + 1] + D[i]
                if i > 0:
                    ksi[i] = C_[i - 1] / (B_[i - 1] - A_[i - 1] * ksi[i - 1])
                    eta[i] = (A_[i - 1] * eta[i - 1] + D_[i - 1]) / (B_[i - 1] - A_[i - 1] * ksi[i - 1])
            for i in range(N, -1, -1):
                if i < N:
                    dy[i] = ksi[i + 1] * dy[i + 1] + eta[i + 1]
                else:
                    dy[i] = (A_[i] * eta[i] + D_[i]) / (B_[i] - A_[i] * ksi[i])
                max_dy = max(max_dy, abs(dy[i] / T[i]))
            if max_dy < EPS:
                break
            else:
                for i in range(N + 1):
                    T[i] += dy[i]
            max_dy = 0
            j += 1

        x = [0 for j in range(N + 1)]
        for i in range(1, N + 1):
            x[i] = x[i - 1] + h
        fig, ax = plt.subplots()
        text = "Зависимость температуры от расстояния от левого торца стержня\n" + "N=" + str(N) + ", F0=" + str(
            F0) + ", a0=" \
               + str(a0) + ", EPS=" + str(EPS)
        plt.title(text)
        ax.plot(x, T, '-')
        ax.set_xlabel('Расстояние от левого торца стержня, см')
        ax.set_ylabel('Температура, К')
        ax.grid(True)
        plt.show()


if __name__ == "__main__":
    main()
