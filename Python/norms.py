import math


def l2_norm(real_values, numerical_values):
    num_values = len(real_values)
    if num_values != len(numerical_values):
        raise ValueError("Real and numerical value lists must be same length")

    combined_values = zip(real_values, numerical_values)
    return math.sqrt(
        1 / num_values
        * sum(
            (real_value - numerical_value) ** 2
            for real_value, numerical_value in combined_values
        ))
