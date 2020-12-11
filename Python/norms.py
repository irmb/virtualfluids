import math


def get_sum_of_squared_distances(real_values, numerical_values):
    combined_values = zip(real_values, numerical_values)
    sum_of_squared_distances = sum((numerical_value - real_value) ** 2
                                   for real_value, numerical_value
                                   in combined_values)
    return sum_of_squared_distances


def root_mean_squared_error(real_values, numerical_values):
    num_values = len(real_values)
    if num_values != len(numerical_values):
        raise ValueError("Real and numerical value lists must be same length")

    sum_of_squared_distances = get_sum_of_squared_distances(real_values, numerical_values)
    sum_of_squared_real_values = sum(real_value ** 2 for real_value in real_values)

    return math.sqrt(sum_of_squared_distances / num_values)


def mean_absolute_error(real_values, numerical_values):
    num_values = len(real_values)
    if num_values != len(numerical_values):
        raise ValueError("Real and numerical value lists must be same length")

    combined_values = zip(real_values, numerical_values)
    sum_of_absolute_distances = sum(abs(numerical_value - real_value)
                                    for real_value, numerical_value
                                    in combined_values)

    return sum_of_absolute_distances / num_values


def mean_squared_error(real_values, numerical_values):
    num_values = len(real_values)
    if num_values != len(numerical_values):
        raise ValueError("Real and numerical value lists must be same length")

    sum_of_squared_distances = get_sum_of_squared_distances(real_values, numerical_values)

    return sum_of_squared_distances / num_values
