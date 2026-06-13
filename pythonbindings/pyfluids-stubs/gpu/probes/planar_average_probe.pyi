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

! \author Henry Korb
=======================================================================================
"""

from __future__ import annotations

import enum
from ... import gpu
from ...basics.geometry3d import Axis

class Statistic(enum.Enum):
    _value_: int
    Means = ...
    Covariances = ...
    Skewness = ...
    Flatness = ...

class PlanarAverageProbe(gpu.Sampler):
    def __init__(
        self,
        para: gpu.Parameter,
        cuda_memory_manager: gpu.CudaMemoryManager,
        probe_name: str,
        output_path: str,
        t_start_avg: int,
        t_start_tmp_avg: int,
        t_avg: int,
        t_start_out: int,
        t_out: int,
        plane_normal: Axis,
        compute_time_averages: bool,
        sample_scalar: bool = ...,
        sample_subgrid_scale_fluxes: bool = ...,
    ) -> None: ...
    def add_statistic(self, variable: Statistic) -> None: ...
    def add_all_available_statistics(self) -> None: ...
    def set_file_name_to_n_out(self) -> None: ...
