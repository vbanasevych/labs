import math

def simpson_method(a, b, epsilon):
    def f(x):
        return 1 / math.sqrt(x**2 + 3.2)

    def simpson(n):
        h = (b - a) / n
        result = f(a) + f(b)

        for i in range(1, n):
            coef = 4 if i % 2 != 0 else 2
            result += coef * f(a + i * h)
        return result * h / 3

    n = 2
    s_n = simpson(n)
    while True:
        n *= 2
        s_2n = simpson(n)
        if abs(s_2n - s_n) < epsilon:
            break
        s_n = s_2n
    return s_2n, n

# Параметри
a = 1.2
b = 2.7
epsilon = 0.01

value, used_n = simpson_method(a, b, epsilon)
print(f"Метод Сімпсона:\nІнтеграл ≈ {value:.6f}, n = {used_n}")
