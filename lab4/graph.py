import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Функція для читання даних з текстового файлу
def read_data(file_path):
    data = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            angle = float(parts[0])   # правильний порядок
            energy = float(parts[1])
            atoms = parts[2]

            if atoms not in data:
                data[atoms] = {'angles': [], 'energies': []}
            data[atoms]['angles'].append(angle)
            data[atoms]['energies'].append(energy)
    return data

# Читання даних з файлу
file_path = r'D:\Аня_лаби\labs\lab4\ehb_results.txt'
data = read_data(file_path)

# Вибір кольорової палітри
num_colors = len(data)
colors = cm.get_cmap('tab20', num_colors)(np.linspace(0, 1, num_colors))

# Створення графіка
plt.figure(figsize=(15, 9))

# Проходимо по кожній групі атомів
for idx, (atoms, values) in enumerate(data.items()):
    plt.scatter(values['angles'], values['energies'], label=atoms, color=colors[idx])

# Додавання підписів та заголовку
plt.title('Залежність енергії водневого зв\'язку від кута обертання')
plt.xlabel('Кут (градуси)')
plt.ylabel('Енергія водневого зв\'язку (ккал/моль)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()

# Показ графіка
plt.show()
