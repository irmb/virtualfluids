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

  You should have received a copy of the GNU General Public License along
  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.

! \file grid_generator.pyi
! \ingroup gpu
! \author Henry Korb
=======================================================================================
"""
from __future__ import annotations

from typing import List

from typing import overload
import gpu



class BoundingBox:
    def __init__(self, min_x: float, max_x: float, min_y: float, max_y: float, min_z: float, max_z: float) -> None: ...


class Object:
    def __init__(self, *args, **kwargs) -> None: ...


class Conglomerate(Object):
    def __init__(self, *args, **kwargs) -> None: ...
    def add(self, object: Object) -> None: ...
    @staticmethod
    def make_shared() -> Conglomerate: ...
    def subtract(self, object: Object) -> None: ...


class Cuboid(Object):
    def __init__(self, min_x1: float, min_x2: float, min_x3: float, max_x1: float, max_x2: float, max_x3: float) -> None: ...


class GridBuilder:
    def __init__(self, *args, **kwargs) -> None: ...
    def get_number_of_grid_levels(self) -> int: ...


class GridFactory:
    def __init__(self, *args, **kwargs) -> None: ...
    @staticmethod
    def make() -> GridFactory: ...


class LevelGridBuilder(GridBuilder):
    def __init__(self, *args, **kwargs) -> None: ...
    def set_no_slip_boundary_condition(self, side_type: gpu.SideType) -> None: ...
    def set_periodic_boundary_condition(self, periodic_x: bool, periodic_y: bool, periodic_z: bool) -> None: ...
    def set_precursor_boundary_condition(self, side_type: gpu.SideType, file_collection: gpu.FileCollection, n_t_read: int, velocity_x: float = ..., velocity_y: float = ..., velocity_z: float = ..., file_level_to_grid_level_map: List[int] = ...) -> None: ...
    def set_pressure_boundary_condition(self, side_type: gpu.SideType, rho: float) -> None: ...
    def set_slip_boundary_condition(self, side_type: gpu.SideType, normal_x: float, normal_y: float, normal_z: float) -> None: ...
    def set_stress_boundary_condition(self, side_type: gpu.SideType, normal_x: float, normal_y: float, normal_z: float, sampling_offset: int, z0: float, dx: float, q: float = ...) -> None: ...
    def set_velocity_boundary_condition(self, side_type: gpu.SideType, vx: float, vy: float, vz: float) -> None: ...



class MultipleGridBuilder(LevelGridBuilder):
    def __init__(self, *args, **kwargs) -> None: ...
    def add_coarse_grid(self, start_x: float, start_y: float, start_z: float, end_x: float, end_y: float, end_z: float, delta: float) -> None: ...
    @overload
    def add_geometry(self, solid_object: Object) -> None: ...
    @overload
    def add_geometry(self, solid_object: Object, level: int) -> None: ...
    @overload
    def add_grid(self, grid_shape: Object) -> None: ...
    @overload
    def add_grid(self, grid_shape: Object, level_fine: int) -> None: ...
    def build_grids(self, enable_thin_walls: bool) -> None: ...
    def get_number_of_levels(self) -> int: ...
    @staticmethod
    def make_shared(grid_factory: GridFactory) -> MultipleGridBuilder: ...


class Sphere(Object):
    def __init__(self, *args, **kwargs) -> None: ...
    @staticmethod
    def make_shared() -> Sphere: ...


class TriangularMesh(Object):
    def __init__(self, *args, **kwargs) -> None: ...
    @staticmethod
    def make() -> TriangularMesh: ...