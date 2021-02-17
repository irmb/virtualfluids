import math
import pyvista as pv


def vertical_column_from_mesh(mesh):
    last_seen = math.inf
    relevant_indices = []
    first_x = 0
    first_y = 0
    for index, point in enumerate(mesh.points):
        if index == 0:
            first_x = point[0]
            first_y = point[1]

        if (point[0] != first_x or point[1] != first_y) and point[2] == last_seen:
            continue

        relevant_indices.append(index)
        last_seen = point[2]

    return relevant_indices


def get_values_from_indices(array, indices):
    return [array[index] for index in indices]


if __name__ == "__main__":
    mesh = pv.read("output/mq/mq10000/mq0_10000.ascii.vtu")
    indices = vertical_column_from_mesh(mesh)
    values = get_values_from_indices(mesh.get_array("Vx"), indices)
    print(len(indices))
    print(values)
