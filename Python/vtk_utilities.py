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
from typing import List

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


def get_values_from_indices(array, indices) -> List[float]:
    return [array[index] for index in indices]


if __name__ == "__main__":
    mesh = pv.read("output/mq/mq10000/mq0_10000.ascii.vtu")
    indices = vertical_column_from_mesh(mesh)
    values = get_values_from_indices(mesh.get_array("Vx"), indices)
    print(len(indices))
    print(values)
