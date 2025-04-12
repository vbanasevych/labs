import numpy as np
import matplotlib.pyplot as plt

# Задання функції f(t, y, y') для рівняння
def f(t, y, y_prime):
    return 4 * y_prime - 13 * y - np.sin(t)

# Ініціалізація початкових умов
t0 = 5
y0 = 7
y_prime0 = 11
h = 0.22  # Крок інтегрування (залежно від кількості кроків)
t_final = 9.5
n_steps = int((t_final - t0) / h)  # Кількість кроків

# Створення масивів для зберігання значень t, y та y'
t_values = np.linspace(t0, t_final, n_steps + 1)
y_values = np.zeros(n_steps + 1)
y_prime_values = np.zeros(n_steps + 1)

# Початкові умови
y_values[0] = y0
y_prime_values[0] = y_prime0

# Основний цикл для методу Хойда
for i in range(n_steps):
    t = t_values[i]
    y = y_values[i]
    y_prime = y_prime_values[i]

    # Оцінка f на поточному кроці
    f_current = f(t, y, y_prime)

    # Прогнозування y' на наступному кроці (метод Ейлера)
    y_prime_predict = y_prime + h * f_current

    # Оцінка f на наступному кроці
    f_next = f(t + h, y + h * y_prime, y_prime + h * f_current)

    # Оновлення значень за методом Хойда
    y_prime_values[i + 1] = y_prime + h * (f_current + f_next) / 2
    y_values[i + 1] = y + h * y_prime_values[i + 1]

# Побудова графіка результатів
plt.plot(t_values, y_values, label="y(t) (метод Хойда)")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Чисельне розв'язання рівняння методом Хойда")
plt.grid(True)
plt.legend()
plt.show()
