# analyze CPS

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def main(cps1_finger: str, cps2_finger: str, directory: str):
	for filename in os.listdir(directory):
		if filename.endswith(".cpsa"):
			print(filename)
			# file = cr39py.CR39(filename)
			# file.add_cut(cr39py.Cut(cmax=max_contrast, emax=max_eccentricity,
			#                         dmin=min_diameter, dmax=max_diameter))
			# x_tracks, y_tracks = file.get_x(), file.get_y()
			# plt.histogram(x_tracks, y_tracks, bins=216)
			# plt.show()

if __name__ == "__main__":
	if len(sys.argv) == 4:
		cps1_finger, cps2_finger, directory = sys.argv[1:]
		with open("arguments.txt", "w") as file:
			file.write(f"cps1-finger={cps1_finger}\ncps2-finger={cps2_finger}\ndirectory={directory}\n")
	elif len(sys.argv) == 1:
		cps1_finger, cps2_finger, directory = None, None, None
		try:
			with open("arguments.txt", "r") as file:
				for line in file:
					if "=" in line:
						key = line[:line.index("=")].strip()
						value = line[line.index("=") + 1:].strip()
						if "cps1" in key:
							cps1_finger = value
						elif "cps2" in key:
							cps2_finger = value
						elif "directory" in key or "path" in key:
							directory = value
						else:
							raise ValueError("The `arguments.txt` file contains an unrecognized key: {key}")
		except IOError:
			raise ValueError("You must run this script with three command line arguments: "
			                 "`python analyze_cps.py cps1-finger cps2-finger directory`; "
			                 "or specify the arguments by creating a file called `arguments.txt` formatted like"
			                 "\n  “cps1-finger=a1\n   cps2-finger=b2w\n   directory=example/path/”")
		if cps1_finger is None or cps2_finger is None or directory is None:
			raise ValueError("The `arguments.txt` file was missing some of the three required arguments.")
	else:
		raise ValueError("You must run this script with exactly three command line arguments: "
		                 "`python analyze_cps.py cps1-finger cps2-finger directory`, "
		                 "or alternatively specify the arguments by creating a file called `arguments.txt` formatted like "
		                 "\n  “cps1-finger=a1\n   cps2-finger=b2w\n   directory=example/path/”")

	main(cps1_finger, cps2_finger, directory)
