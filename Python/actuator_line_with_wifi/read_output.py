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
#%%
import tables as tb
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from wiFI.turbine import create_turbine_from_farm_json


def average_data(file: tb.File, key: str):
    return np.array(file.root[key]["averages"].data)

def instantaneous_data(file: tb.File, key: str):
    return np.array(file.root[key]["instantaneous"].data)

radius = 63
n_nodes = 32
velocity = 9
#%%
turbine_file = Path(__file__).parent/"SingleTurbine.json"
output_dir = Path(__file__).parent/"output"
wifi_output = tb.File(output_dir/"wifi.h5", mode="r")
turbine_model = create_turbine_from_farm_json(turbine_file, number_of_blade_nodes=32)
# %%
print(wifi_output.root)
#%%
rotor_speed = pd.DataFrame(wifi_output.root.rotor_speed.averages.data, wifi_output.root.rotor_speed.averages.times, copy=False)
plt.plot(rotor_speed)
# %%
blade_forces = average_data(wifi_output, "blade_forces")
fig, ax = plt.subplots(3, sharex=True)
force_factor = 1.225*radius*velocity**2
dr = radius/n_nodes
ax[0].plot(blade_forces[-1, 0, 0, :, :].T/dr/force_factor)
ax[1].plot(-blade_forces[-1, 1, 0, :, :].T/dr/force_factor)
ax[2].plot(-blade_forces[-1, 2, 0, :, :].T/dr/force_factor)

# %%
blade_coordinates = instantaneous_data(wifi_output, "blade_coordinates")
for i in range(3):
    plt.scatter(blade_coordinates[-1,1,0,i,:], blade_coordinates[-1,2,0,i,:])

# %%
blade_velocities = average_data(wifi_output, "blade_velocities")
plt.plot(blade_velocities[-1, 0, 0, :, :].T)
# %%
fig, ax = plt.subplots(3, sharex=True)
for i, axis in enumerate(ax):
    axis.plot(blade_velocities[-1, i, 0,:,:].T)

# %%
