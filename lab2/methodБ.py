import math

def f(x):
    return 1 / math.sqrt(x ** 2 + 3.2)

def right_rectangle_method(a, b, epsilon):
    n = 2
    prev_result = 0

    while True:
        h = (b - a) / n
        result = 0
        for i in range(1, n + 1):
            xi = a + i * h
            result += f(xi)
        result *= h

        if abs(result - prev_result) < epsilon:
            break

        prev_result = result
        n *= 2
    return result, n

# Параметри
a = 1.2
b = 2.7
epsilon = 0.001

value, used_n = right_rectangle_method(a, b, epsilon)
print(f"Метод правих прямокутників:\nІнтеграл ≈ {value:.6f}, n = {used_n}")