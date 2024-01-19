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
from dataclasses import dataclass


@dataclass
class PoiseuilleSettings:
    density = 1
    viscosity = 0.005
    height = 10
    length = 1
    pressure_in = 0
    pressure_out = 0
    force = 0


def poiseuille_at_z(settings: PoiseuilleSettings, z: float):
    pressure_grad = ((settings.pressure_out - settings.pressure_in) / settings.length)

    return ((1 / settings.viscosity)
            * (- pressure_grad + settings.density * settings.force)
            * (z / 2)
            * (settings.height - z))


def poiseuille_at_heights(settings: PoiseuilleSettings, heights):
    return [poiseuille_at_z(settings, z) for z in heights]


def reynolds_number(settings: PoiseuilleSettings):
    max_v = poiseuille_at_z(settings, settings.height / 2)
    return max_v * settings.height / settings.viscosity


if __name__ == '__main__':
    sim_settings = PoiseuilleSettings()

    sim_settings.force = 2e-7
    sim_settings.viscosity = 1e-3
    sim_settings.height = 16
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))

    sim_settings.viscosity *= 2
    sim_settings.height *= 2
    sim_settings.force /= 2
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))

    sim_settings.viscosity *= 2
    sim_settings.height *= 2
    sim_settings.force /= 2
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))
