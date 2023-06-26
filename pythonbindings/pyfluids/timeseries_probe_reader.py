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