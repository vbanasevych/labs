import math

def f(x):
    return 3 * x - math.cos(x) - 1

def bisection_method(f, a, b, eps):
    iterations = 0
    x_prev2 = a
    x_prev1 = b
    while (b - a) / 2 > eps:
        c = (a + b) / 2
        fc = f(c)
        iterations += 1

        # Умова (9.16)
        if abs(fc) < eps:
            return c, iterations

        # Умова (9.18) — ітераційне уточнення
        if iterations >= 3:
            num = (c - x_prev1) ** 2
            denom = 2 * abs(c - x_prev1 - x_prev2)
            if denom != 0 and num / denom < eps:
                return c, iterations

        if f(a) * fc < 0:
            b = c
        else:
            a = c

        x_prev2 = x_prev1
        x_prev1 = c

    return (a + b) / 2, iterations

def chord_method(f, a, b, eps):
    x_prev2 = a
    x_prev1 = b
    iterations = 0
    while True:
        fx = f(x_prev1)
        fx_prev = f(x_prev2)
        x_new = x_prev1 - fx * (x_prev1 - x_prev2) / (fx - fx_prev)
        iterations += 1

        # Умова (9.16)
        if abs(f(x_new)) < eps:
            return x_new, iterations

        # Умова (9.18)
        if iterations >= 3:
            num = (x_new - x_prev1) ** 2
            denom = 2 * abs(x_new - x_prev1 - x_prev2)
            if denom != 0 and num / denom < eps:
                return x_new, iterations

        x_prev2 = x_prev1
        x_prev1 = x_new

# Початкові параметри
a, b = 0, 1
eps = 0.01

# Виклик методів
root_bisection, iter_bis = bisection_method(f, a, b, eps)
root_chord, iter_chord = chord_method(f, a, b, eps)

# Вивід результатів
print("Метод половинного ділення:")
print(f"Корінь: {root_bisection:.6f}, f(x): {f(root_bisection):.6f}, ітерацій: {iter_bis}")

print("\nМетод хорд:")
print(f"Корінь: {root_chord:.6f}, f(x): {f(root_chord):.6f}, ітерацій: {iter_chord}")
