import numpy as np
import os

periodic_table = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
    19: 'K', 20: 'Ca'
}

def read_xyz(file_path):
    """Зчитування xyz-файлу та отримання списку атомів і координат."""
    atoms = []
    coordinates = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines:
        parts = line.split()
        if len(parts) == 4:
            atom_number = int(parts[0])
            atom = periodic_table.get(atom_number, str(atom_number))
            x, y, z = map(float, parts[1:])
            atoms.append(atom)
            coordinates.append([x, y, z])
    return atoms, np.array(coordinates)


def calculate_distance(coord1, coord2):
    """Обчислення відстані між двома атомами."""
    return np.linalg.norm(coord1 - coord2)


def build_bond_matrix(atoms, coordinates, threshold=1.6):
    """Побудова матриці зв'язків на основі відстаней між атомами.
    Гідроген не може утворювати зв'язок з іншим гідрогеном."""
    num_atoms = len(atoms)
    bond_matrix = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            # Пропускаємо випадок, коли обидва атоми - гідроген
            if atoms[i] == 'H' and atoms[j] == 'H':
                continue

            # Обчислюємо відстань між атомами
            distance = calculate_distance(coordinates[i], coordinates[j])

            # Якщо відстань менша за порогове значення, створюємо зв'язок
            if distance < threshold:
                bond_matrix[i, j] = 1
                bond_matrix[j, i] = 1

    return bond_matrix


def find_cycles(bond_matrix):
    """Пошук циклів у молекулі за допомогою пошуку в глибину."""
    num_atoms = len(bond_matrix)
    visited = [False] * num_atoms
    cycles = []

    def dfs(current, parent, path):
        visited[current] = True
        path.append(current)

        for neighbor in np.nonzero(bond_matrix[current])[0]:
            if neighbor == parent:
                continue
            if neighbor in path:
                cycle = path[path.index(neighbor):]
                if cycle not in cycles:
                    cycles.append(cycle)
            elif not visited[neighbor]:
                dfs(neighbor, current, path)

        path.pop()

    for i in range(num_atoms):
        if not visited[i]:
            dfs(i, -1, [])

    return cycles


def find_glycosidic_bond(bond_matrix, atoms, cycles):
    """Пошук глікозидного зв'язку O-C-N-C з такими умовами:
    1. Атом кисню і перший карбон мають бути в одному циклі.
    2. Атом нітрогену має бути в іншому циклі, ніж перші два атоми, але в тому ж циклі, що й останній карбон.
    3. Останній карбон не повинен мати сусідів-гідрогенів.
    """

    # Проходимо по всіх циклах для пошуку атома кисню
    for cycle1 in cycles:
        for atom1 in cycle1:
            if atoms[atom1] == 'O':

                # Знаходимо сусідній атом карбону (atom2)
                neighbors_c1 = np.nonzero(bond_matrix[atom1])[0]
                for atom2 in neighbors_c1:
                    if atoms[atom2] == 'C' and atom2 in cycle1:  # Карбон у тому ж циклі, що й кисень

                        # Шукаємо сусідній атом нітрогену (atom3)
                        neighbors_c2 = np.nonzero(bond_matrix[atom2])[0]
                        for atom3 in neighbors_c2:
                            if atoms[atom3] == 'N':
                                # Нітроген має бути у іншому циклі, ніж O і C
                                cycle2 = next((cycle for cycle in cycles if atom3 in cycle and cycle != cycle1), None)
                                if not cycle2:
                                    continue

                                # Шукаємо останній атом карбону (atom4), який зв'язаний із нітрогеном (atom3)
                                neighbors_c3 = np.nonzero(bond_matrix[atom3])[0]
                                for atom4 in neighbors_c3:
                                    if atoms[atom4] == 'C' and atom4 in cycle2:
                                        # Перевіряємо, що цей атом карбону не має сусідів-гідрогенів
                                        neighbors_c4 = np.nonzero(bond_matrix[atom4])[0]
                                        if all(atoms[neighbor] != 'H' for neighbor in neighbors_c4):

                                            # Переконуємося, що атоми пов'язані ланцюжком O-C-N-C
                                            if (bond_matrix[atom1, atom2] and
                                                    bond_matrix[atom2, atom3] and
                                                    bond_matrix[atom3, atom4]):
                                                print(
                                                    f"Знайдений глікозидний кут O-C-N-C: {atoms[atom1]}-{atoms[atom2]}-{atoms[atom3]}-{atoms[atom4]}")
                                                print(
                                                    f"Індекси атомів: {atom1 + 1}, {atom2 + 1}, {atom3 + 1}, {atom4 + 1}")
                                                return atom1, atom2, atom3, atom4
    return None


def find_atoms_to_rotate(coordinates, c2, c3):
    """
    Знаходження атомів для обертання на основі відстані до атома з більшим індексом.
    """
    # Визначаємо ближній і дальній атоми за індексом
    closer_atom, farther_atom = (c2, c3) if c2 > c3 else (c3, c2)

    atoms_to_rotate = []

    # Обчислюємо відстані до ближнього та дальнього атомів для кожного атома
    for atom_id in range(len(coordinates)):
        if atom_id in (c2, c3):
            continue  # Пропускаємо самі атоми глікозидного зв'язку

        # Розраховуємо відстані до ближнього і дальнього атомів
        distance_to_closer = np.linalg.norm(coordinates[atom_id] - coordinates[closer_atom])
        distance_to_farther = np.linalg.norm(coordinates[atom_id] - coordinates[farther_atom])

        # Додаємо атом до списку, якщо він ближчий до ближнього атома
        if distance_to_closer < distance_to_farther:
            atoms_to_rotate.append(atom_id)

    return atoms_to_rotate


def rotate_coordinates(coords, atoms_to_rotate, point, axis, angle):
    """Обертання координат вибраних атомів навколо осі на заданий кут."""
    angle_rad = np.radians(angle)
    axis = axis / np.linalg.norm(axis)

    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)

    ux, uy, uz = axis
    rotation_matrix = np.array([
        [cos_angle + ux ** 2 * (1 - cos_angle), ux * uy * (1 - cos_angle) - uz * sin_angle,
         ux * uz * (1 - cos_angle) + uy * sin_angle],
        [uy * ux * (1 - cos_angle) + uz * sin_angle, cos_angle + uy ** 2 * (1 - cos_angle),
         uy * uz * (1 - cos_angle) - ux * sin_angle],
        [uz * ux * (1 - cos_angle) - uy * sin_angle, uz * uy * (1 - cos_angle) + ux * sin_angle,
         cos_angle + uz ** 2 * (1 - cos_angle)]
    ])

    rotated_coords = coords.copy()
    for atom_index in atoms_to_rotate:
        shifted_coord = coords[atom_index] - point
        rotated_coord = np.dot(rotation_matrix, shifted_coord) + point
        rotated_coords[atom_index] = rotated_coord
    return rotated_coords


def generate_gaussian_files(atoms, coordinates, torsion_atoms, atoms_to_rotate):
    """Генерація вихідних файлів для обертання на різні кути."""
    output_dir = "gaussian_files"
    os.makedirs(output_dir, exist_ok=True)

    for angle in range(0, 360):
        rotated_coords = rotate_coordinates(
            coordinates, atoms_to_rotate,
            coordinates[torsion_atoms[1]],
            coordinates[torsion_atoms[2]] - coordinates[torsion_atoms[1]],
            angle
        )
        file_name = os.path.join(output_dir, f"mol_{angle}.gjf")
        with open(file_name, 'w') as file:
            # file.write(f"{len(atoms)}\n\n")
            file.write("# hf/3-21+G Opt=ModRedundant output=wfn\n\n")
            file.write(f"Torsion angle {angle} degrees\n\n0 1\n")
            for atom, (x, y, z) in zip(atoms, rotated_coords):
                file.write(f"{atom:2} {x:15.9f} {y:15.9f} {z:15.9f}\n")
            file.write("\n")
            file.write(f"D {c1 + 1} {c2 + 1} {c3 + 1} {c4 + 1} F\n\n")
            file.write(f"result_{angle}.wfn")
        print(f"Створено файл: {file_name}")


# Основна частина
input_file = "D:\Аня_лаби\mol.xyz"
atoms, coordinates = read_xyz(input_file)

bond_matrix = build_bond_matrix(atoms, coordinates)
cycles = find_cycles(bond_matrix)

glycosidic_bond = find_glycosidic_bond(bond_matrix, atoms, cycles)

if glycosidic_bond:
    c1, c2, c3, c4 = glycosidic_bond
    print(f"Глікозидний зв'язок знайдено між атомами: {atoms[c2]}-{atoms[c3]} ({c2 + 1}, {c3 + 1})")

    atoms_to_rotate = find_atoms_to_rotate(coordinates, c2, c3)
    if atoms_to_rotate:
        generate_gaussian_files(atoms, coordinates, (c1, c2, c3, c4), atoms_to_rotate)
    else:
        print("Не вдалося знайти атоми для обертання.")
else:
    print("Глікозидний зв'язок між циклами не знайдено.")
