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
