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
import numpy as np
import scipy.stats as stats
import errors
from SlurmTests.poiseuille.result_collector import collect_results
from SlurmTests.poiseuille.settings import Scaling

analytical_results, numerical_results = collect_results()
normalized_l2_errors = [errors.normalized_l2_error(analytical, numerical)
                        for analytical, numerical in zip(analytical_results, numerical_results)]

nodes_in_x3_per_run = []
for simulation_run in range(0, 3):
    grid_params, _, _, _ = Scaling.configuration_for_scale_level(simulation_run)
    nodes_in_x3_per_run.append(grid_params.number_of_nodes_per_direction[2])

nodes_as_log = [np.log10(node) for node in nodes_in_x3_per_run]
l2_norms_as_log = [np.log10(l2) for l2 in normalized_l2_errors]
res = stats.linregress(nodes_as_log, l2_norms_as_log)

assert res.slope <= -2, f"Expected slope of l2 error to be <= -2, but was {res.slope}"
