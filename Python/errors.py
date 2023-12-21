r"""
=======================================================================================
 ____          ____    __    ______     __________   __      __       __        __
 \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
      \    \  |    |   ________________________________________________________________
       \    \ |    |  |  ______________________________________________________________|
        \    \|    |  |  |         __          __     __     __     ______      _______
         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/

  This file is part of VirtualFluids. VirtualFluids is free software: you can
  redistribute it and/or modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  SPDX-License-Identifier: GPL-3.0-or-later
  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

=======================================================================================
"""
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


def normalized_l2_error(real_values, numerical_values):
    sum_of_squared_distances = get_sum_of_squared_distances(real_values, numerical_values)
    sum_of_squared_real_values = sum(real_value ** 2 for real_value in real_values)

    return math.sqrt(sum_of_squared_distances / sum_of_squared_real_values)
