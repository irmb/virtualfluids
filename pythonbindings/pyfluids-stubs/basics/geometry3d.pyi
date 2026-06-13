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
from enum import Enum
from numpy.typing import NDArray
import numpy as np

class Axis(Enum):
    x = ...
    y = ...
    z = ...

class GbSpatialData3DReal: ...

class GbStructuredMesh3DReal(GbSpatialData3DReal):
    class InterpolationStrategy(Enum):
        Trilinear = ...
    def __init__(
        self,
        spacing: NDArray[np.float32],
        origin: NDArray[np.float32],
        n_points: NDArray[np.uint32],
        values: NDArray[np.float32],
        interpolation_strategy: GbStructuredMesh3DReal.InterpolationStrategy,
    ): ...

class GbPointCloud3DReal(GbSpatialData3DReal):
    class InterpolationStrategy: ...
    class InverseDistanceWeighingReal(GbPointCloud3DReal.InterpolationStrategy):
        @staticmethod
        def make(n_points: int) -> GbPointCloud3DReal.InterpolationStrategy: ...
    class NearestNeighborReal(GbPointCloud3DReal.InterpolationStrategy):
        @staticmethod
        def make() -> GbPointCloud3DReal.InterpolationStrategy: ...

    def __init__(
        self,
        coords: NDArray[np.float32],
        data: NDArray[np.float32],
        interpolation_strategy: GbPointCloud3DReal.InterpolationStrategy,
        print_tree: bool = ...,
    ): ...

class GbSpatialData3DReal3: ...

class GbStructuredMesh3DReal3(GbSpatialData3DReal3):
    class InterpolationStrategy(Enum):
        Trilinear = ...
    def __init__(
        self,
        spacing: NDArray[np.float32],
        origin: NDArray[np.float32],
        n_points: NDArray[np.uint32],
        values: NDArray[np.float32],
        interpolation_strategy: GbStructuredMesh3DReal3.InterpolationStrategy,
    ): ...

class GbPointCloud3DReal3(GbSpatialData3DReal3):
    class InterpolationStrategy: ...
    class InverseDistanceWeighingReal(GbPointCloud3DReal3.InterpolationStrategy):
        @staticmethod
        def make(n_points: int) -> GbPointCloud3DReal3.InterpolationStrategy: ...
    class NearestNeighborReal(GbPointCloud3DReal3.InterpolationStrategy):
        @staticmethod
        def make() -> GbPointCloud3DReal3.InterpolationStrategy: ...

    def __init__(
        self,
        coords: NDArray[np.float32],
        data: NDArray[np.float32],
        interpolation_strategy: GbPointCloud3DReal3.InterpolationStrategy,
        print_tree: bool = ...,
    ): ...