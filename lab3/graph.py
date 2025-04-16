import numpy as np
import matplotlib.pyplot as plt

# Функція
def f(x):
    return 3 * x - np.cos(x) - 1

# Для візуалізації
x_vals = np.linspace(0, 1.5, 400)
y_vals = f(x_vals)

# Побудова графіка
plt.figure(figsize=(8, 5))
plt.plot(x_vals, y_vals, label=r"$f(x) = 3x - \cos{x} - 1$")
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.title("Графік функції $f(x) = 3x - \cos{x} - 1$")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()
plt.show()
