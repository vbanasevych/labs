import math

# Функція
def f(x):
    return 3*x - math.cos(x) - 1

# Метод А: Половинного ділення з протоколом
def bisection_method(a, b, eps):
    iterations = 0
    history = []
    while (b - a)/2 > eps:
        c = (a + b)/2
        fc = f(c)
        history.append((iterations + 1, a, b, c, fc))
        if abs(fc) < eps:
            break
        if f(a)*fc < 0:
            b = c
        else:
            a = c
        iterations += 1
    return (a + b)/2, iterations, history

# Метод Б2: Хорд з протоколом
def chord_method(a, b, eps):
    x_prev = a
    x = b
    iterations = 0
    history = [(iterations + 1, x_prev, x, f(x))]
    while abs(x - x_prev) > eps or abs(f(x)) > eps:
        x_new = x - f(x)*(x - x_prev)/(f(x) - f(x_prev))
        x_prev, x = x, x_new
        iterations += 1
        history.append((iterations + 1, x_prev, x, f(x)))
    return x, iterations, history

# Початкові умови
a, b = 0.6, 0.7
eps = 0.01

# Обчислення
bisect_root, bisect_iter, bisect_table = bisection_method(a, b, eps)
chord_root, chord_iter, chord_table = chord_method(a, b, eps)

# Вивід результатів
print("Метод половинного ділення (А):")
print(f"Корінь: {bisect_root:.5f}, Ітерацій: {bisect_iter}")
print("Таблиця ітерацій:")
print("№\t\t a\t\t\t b\t\t\t c\t\t f(c)")
for i, a_, b_, c_, fc in bisect_table:
    print(f"{i}\t {a_:.5f}\t {b_:.5f}\t {c_:.5f}\t {fc:.5f}")

print("\nМетод хорд (Б2):")
print(f"Корінь: {chord_root:.5f}, Ітерацій: {chord_iter}")
print("Таблиця ітерацій:")
print("№\t x_prev\t\t x\t\t f(x)")
for i, x_prev, x, fx in chord_table:
    print(f"{i}\t {x_prev:.5f}\t {x:.5f}\t {fx:.5f}")
