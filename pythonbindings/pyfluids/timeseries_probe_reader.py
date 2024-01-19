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
import numpy as np
from pathlib import Path
import pandas as pd
#%%


class TimeseriesProbeReader:
    def __init__(self, file: Path):
        self.file = file
        self.quants, self.positions, self.data = \
            self.read_file()

    def read_file(self):
        with open(self.file, "rb") as f:
            header_length = 0
            header_length += len(f.readline()) # first line
            quant_line = f.readline()
            header_length += len(f.readline()) # number of points
            number_of_points_line = f.readline()
            header_length += len(f.readline()) # positions
            n_points = int(number_of_points_line)
            positions = np.zeros((n_points, 3))
            for i in range(n_points):
                pos_line = f.readline()
                header_length += len(pos_line)
                positions[i] = [float(pos) for pos in pos_line.split(b", ")]

        header_length += len(quant_line)
        header_length += len(number_of_points_line)

        quants = quant_line.decode().split(" ")[1:-1]
        n_quants = len(quants)
        data = np.fromfile(self.file, dtype=np.float32, offset=header_length)
        n_timesteps = len(data)//(n_quants*n_points+1)
        return quants, positions, data.reshape(n_timesteps, n_points*n_quants+1)
    
    def get_data(self):
        return self.data
    
    def get_positions(self):
        return self.positions
    
    def get_quantities(self):
        return self.quants
    
    def to_dataframe(self):
        return pd.DataFrame(self.data[:,1:], columns=self.quants, index=self.data[:,0])