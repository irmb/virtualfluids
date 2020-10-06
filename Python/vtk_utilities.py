import math
import pyvista as pv


def vertical_column_from_mesh(mesh):
    last_seen = math.inf
    relevant_indices = []
    for index, point in enumerate(mesh.points):
        if point[2] == last_seen:
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
