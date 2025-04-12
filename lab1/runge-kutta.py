import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Початкові умови та інтервал
t0, t1 = 5, 9.5
y0, dy0 = 7, 11
steps = 20
t_eval = np.linspace(t0, t1, steps + 1)

# Система ОДУ (переписане з другого порядку у перший)
def system(t, Y):
    y, dy = Y
    ddy = 4*dy - 13*y - np.sin(t)
    return [dy, ddy]

# Чисельне розв'язання
sol = solve_ivp(system, [t0, t1], [y0, dy0], t_eval=t_eval)

# Побудова графіка
plt.figure(figsize=(10, 5))
plt.plot(sol.t, sol.y[0], 'o-', color='blue', label='y(t)')
plt.title("Розв'язок рівняння y'' - 4y' + 13y = -sin(t)")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.grid(True)
plt.legend()
plt.show()

# Таблиця значень y(t)
print(f"{'t':>10} | {'y(t)':>15}")
print("-" * 28)
for t, y_val in zip(sol.t, sol.y[0]):
    print(f"{t:10.3f} | {y_val:15.6f}")
