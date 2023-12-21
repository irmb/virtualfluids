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
from pandas import HDFStore
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

from wiFI.turbine import create_turbine_from_farm_json

def get_blade_array(hdf5: HDFStore, key: str):
    blade_df = hdf5.get(key)
    _, n_turbines, n_blades, n_nodes = blade_df.columns[-1]
    return blade_df.to_numpy().reshape(-1, 3, n_turbines+1, n_blades+1, n_nodes+1)

#%%
turbine_file = Path(__file__).parent/"SingleTurbine.json"
output_dir = Path(__file__).parents[2]/"output"
wifi_output = HDFStore(output_dir/"wifi.h5", mode="r")
turbine_model = create_turbine_from_farm_json(turbine_file, number_of_blade_nodes=32)
# %%
print(wifi_output.keys())
#%%

plt.plot(wifi_output["/Rotor_speed"])
# %%
blade_forces = get_blade_array(wifi_output, "Blade_forces")

plt.plot(blade_forces[-1, 0, 0, :, :].T)

# %%
blade_coordinates = get_blade_array(wifi_output, "Blade_coordinates")
plt.scatter(blade_coordinates[-1,1,0,:,:], blade_coordinates[-1,2,0,:,:])

# %%
blade_velocities = get_blade_array(wifi_output, "Blade_velocities")
plt.plot(blade_velocities[-1, 0, 0, :, :].T)
# %%
plt.plot(blade_velocities[-1, 1, 0,:,:].T)

# %%
