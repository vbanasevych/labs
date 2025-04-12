import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Початкові умови та інтервал
t0, t1 = 5, 9.5
y0, dy0 = 7, 11
steps = 20
t_eval = np.linspace(t0, t1, steps + 1)
h = (t1 - t0) / steps  # Крок для методу Ейлера

# Система ОДУ (переписане з другого порядку у перший)
def system(t, Y):
    y, dy = Y
    ddy = 4*dy - 13*y - np.sin(t)
    return [dy, ddy]

# Чисельне розв'язання методом Рунге-Кутти
sol_rk = solve_ivp(system, [t0, t1], [y0, dy0], t_eval=t_eval)

# Метод Ейлера
def euler_method(f, t0, y0, dy0, h, steps):
    t_values = np.linspace(t0, t1, steps + 1)
    y_values = np.zeros(steps + 1)
    dy_values = np.zeros(steps + 1)
    y_values[0], dy_values[0] = y0, dy0

    for i in range(steps):
        t_values[i + 1] = t_values[i] + h
        dy_values[i + 1] = dy_values[i] + h * (4 * dy_values[i] - 13 * y_values[i] - np.sin(t_values[i]))
        y_values[i + 1] = y_values[i] + h * dy_values[i]

    return t_values, y_values

# Розв'язок методом Ейлера
t_euler, y_euler = euler_method(system, t0, y0, dy0, h, steps)

# Побудова графіка
plt.figure(figsize=(10, 5))
plt.plot(sol_rk.t, sol_rk.y[0], 'o-', label="Метод Рунге-Кутти", color='blue')
plt.plot(t_euler, y_euler, 's--', label="Метод Ейлера", color='red')
plt.title("Порівняння методів: Рунге-Кутти та Ейлера")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True)
plt.show()

# Похибка між методами
error = np.abs(sol_rk.y[0] - np.interp(sol_rk.t, t_euler, y_euler))

# Графік похибки
plt.figure(figsize=(10, 5))
plt.plot(sol_rk.t, error, 'o-', color='green', label='Похибка між методами')
plt.title("Графік похибки між методами Рунге-Кутти та Ейлера")
plt.xlabel("t")
plt.ylabel("Похибка |RK - Euler|")
plt.grid(True)
plt.legend()
plt.show()

# Таблиця значень для методу Рунге-Кутти та Ейлера
print(f"{'t':>10} | {'y(RK)':>15} | {'y(Euler)':>15} | {'Error':>15}")
print("-" * 60)
for t, y_rk, y_euler, err in zip(sol_rk.t, sol_rk.y[0], y_euler, error):
    print(f"{t:10.3f} | {y_rk:15.6f} | {y_euler:15.6f} | {err:15.6f}")
