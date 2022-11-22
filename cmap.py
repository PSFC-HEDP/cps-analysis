import os

import numpy as np
from matplotlib.colors import ListedColormap

CMAP: dict[str, ListedColormap] = {}

for filename in os.listdir("."):
	if filename.startswith("cmap_") and filename.endswith(".csv"):
		cmap_data = np.loadtxt(filename, delimiter=",")
		cmap_name = filename[5:-4]
		CMAP[cmap_name] = ListedColormap(cmap_data, name=cmap_name)
