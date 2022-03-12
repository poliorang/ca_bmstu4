from funcs import *

if __name__ == "__main__":
    # filename = input("Enter filename: ")
    # plist = read_from_file(filename)
    data = read_from_file("../data/data_02.txt")

    if not data:
        print("Ошибка загрузки данных из файла")
        exit(1)

    data.sort(key=lambda x: x[0])

    methods_comparison(data)
    # все исходные игреки, потому что коэффициент a[i] == y[i-1]
    a = [0] + [p[1] for p in data]

    # разницы y[i] - y[i-1]
    h = [0] + [data[i][0] - data[i - 1][0] for i in range(1, len(data))]

    c = [0 for i in range(len(a))]
    b = c[:]
    d = c[:]

    # прогоночные коэффициенты для метода гаусса
    ksi = c[:]
    eta = c[:]

    print("\nРешение задачи")
    x = input_x()
    i = find_section(data, x)
    print("x = ", x, "  ->  x ∈ [", data[i - 1][0], '; ', data[i][0], ']\n', sep='')

    c, ksi, eta = start_ratios(x, data, c, ksi, eta)
    # print(c, ksi, eta, sep='\n')

    straight_walk(ksi, eta, len(c), data, h)
    forward_walk(ksi, eta, c)
    find_ratios(d, b, a, c, h)

    result = find_result(x, a[i], b[i], c[i], d[i], data[i - 1][0])
    print("\nResult -->", result)
