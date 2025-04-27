import os
import re

# Функція для обробки одного файлу
def process_file(file_name):
    with open(file_name, 'r') as file:
        content = file.read()

    # Шукаємо блоки між "Point is a Bond Critical Point (BCP)" і "DelSq(K(r))"
    blocks = re.findall(r"Point is a Bond Critical Point \(BCP\)(.*?)DelSq\(K\(r\)\)", content, re.DOTALL)

    # Списки для збереження результатів
    ehb_results = []

    # Обробка кожного блоку
    for i, block in enumerate(blocks):
        if "Bond path linked to nuclear attractor H" in block:
            # Шукаємо числове значення після "DelSq(Rho(r))" (з урахуванням експоненційної частини)
            rho_match = re.search(r"DelSq\(Rho\(r\)\)\s*([+-]?\d*\.\d+([eE][+-]?\d+)?)", block)
            if rho_match:
                rho_value = float(rho_match.group(1))
                if rho_value > 0:  # тільки якщо позитивне значення
                    # Шукаємо залишки рядків з "Bond path linked to nuclear attractor"
                    atom_matches = re.findall(r"Bond path linked to nuclear attractor\s*(\S+)", block)
                    atom_1 = atom_matches[0] if len(atom_matches) > 0 else "None"
                    atom_2 = atom_matches[1] if len(atom_matches) > 1 else "None"
                    # Шукаємо числове значення після "V(r)" (з урахуванням експоненційної частини)
                    v_match = re.search(r"V\(r\)\s*([+-]?\d*\.\d+([eE][+-]?\d+)?)", block)
                    if v_match:
                        v_value = float(v_match.group(1))
                        v_kcal = v_value * 627.509
                        ehb = 0.5 * v_kcal
                        # Витягуємо кут з назви файлу
                        angle = int(re.search(r'(\d+)', os.path.basename(file_name)).group(1))
                        ehb_results.append(f"{angle} {ehb:.6f} {atom_1}—{atom_2}")

    return ehb_results


# Функція для обробки всіх файлів
def process_all_files():
    # Задаємо шлях до папки з файлами
    folder_path = r'D:\Аня_лаби\labs\lab4\extout_files'

    # Збираємо всі файли .extout в директорії
    files = [f for f in os.listdir(folder_path) if f.endswith('.extout')]

    # Сортуємо файли за числовим значенням кута
    files.sort(key=lambda f: int(re.search(r'(\d+)', f).group(1)))

    # Файл для збереження результатів
    output_path = r'D:\Аня_лаби\labs\lab4\ehb_results.txt'

    with open(output_path, 'w') as ehb_file:
        for file in files:
            file_path = os.path.join(folder_path, file)
            ehb_results = process_file(file_path)

            for result in ehb_results:
                ehb_file.write(result + '\n')
    print("Файл створено.")


# Запуск
process_all_files()