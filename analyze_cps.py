# analyze CPS

import os
import re
import sys
from math import log, floor, ceil

from cr39py import cr39
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from numpy.typing import NDArray

from cmap import CMAP

SLIT_WIDTH = .2  # (cm)
MAX_CONTRAST = 35
MAX_ECCENTRICITY = 15
MAX_DIAMETER = 30

def main(cps1_finger: str, cps2_finger: str, directory: str) -> None:
	for filename in os.listdir(directory):
		if filename.endswith(".cpsa"):
			calibration_x, calibration_min_energy, calibration_nominal_energy, calibration_max_energy = load_calibration(filename, cps1_finger, cps2_finger)
			tracks_x, tracks_y, tracks_d, tracks_c = load_tracks(directory, filename)
			plt.figure()
			histogram2d("x (cm)", tracks_x, "y (cm)", tracks_y, filename[:-4])
			plt.figure()
			histogram2d("x (cm)", tracks_x, "d (μm)", tracks_d, filename[:-4])
			plt.figure()
			histogram2d("d (μm)", tracks_d, "c (%)", tracks_c, filename[:-4], log_scale=True)
			plt.figure()
			plt.fill_between(calibration_x, calibration_min_energy, calibration_max_energy, alpha=.5)
			plt.plot(calibration_x, calibration_nominal_energy)
			plt.xlabel("x (cm)")
			plt.ylabel("Proton energy (MeV)")
			plt.show()


def load_calibration(filename: str, cps1_finger: str, cps2_finger: str) -> tuple[NDArray[float], NDArray[float], NDArray[float]]:
	if "cps1" in filename.lower():
		cps = 1
	elif "cps2" in filename.lower():
		cps = 2
	elif cps2_finger.lower().startswith("n"):
		cps = 1
	elif cps1_finger.lower().startswith("n"):
		cps = 2
	else:
		raise ValueError(f"the filename doesn't make it clear whether CPS1 or CPS2 was used: `{filename}`")
	finger = cps1_finger if cps == 1 else cps2_finger
	x, energy = None, None
	i, j = 0, 0
	try:
		with open(os.path.join("calibrations", f"cps{cps}-{finger}.csv"), "r") as file:
			for line in file:
				if re.fullmatch(r"\d+\s*", line) and x is None:
					resolution = (int(line) - 2)//3
					x = np.empty(resolution)
					energy = np.empty((3, resolution))
				elif re.fullmatch(r"[-+.\de]+ , [-+.\de]+\s*", line):
					values = [float(token) for token in line.split(",")]
					x[j] = values[0]
					energy[i, j] = values[1]
					j += 1
				elif re.fullmatch(r"1 ,\s*", line):
					i += 1
					j = 0
	except IOError:
		raise IOError(f"I couldn't find the calibration file for finger {finger} on CPS{cps}.  "
		              f"please get the calibration from Fredrick’s AnalyzeCR39 program and save it "
		              f"to a file in the `calibrations` directory called `cps{cps}-{finger}.csv`.")
	if x is None or energy is None:
		raise RuntimeError
	return x, energy[0, :], energy[1, :], energy[2, :]


def load_tracks(directory: str, filename: str) -> tuple[NDArray[float], NDArray[float], NDArray[float], NDArray[float]]:
	file = cr39.CR39(os.path.join(directory, filename))
	file.add_cut(cr39.Cut(cmin=MAX_CONTRAST))
	file.add_cut(cr39.Cut(emin=MAX_ECCENTRICITY))
	file.add_cut(cr39.Cut(dmin=MAX_DIAMETER))
	file.add_cut(cr39.Cut(ymin=1.35))
	file.apply_cuts()
	if file.ntracks == 0:
		raise ValueError("the file is now empty")
	tracks_x, tracks_y, tracks_d, tracks_c = file.trackdata_subset[:, [0, 1, 2, 3]].T
	return tracks_x, tracks_y, tracks_d, tracks_c


def histogram2d(x_label: str, x_values: NDArray[float], y_label: str, y_values: NDArray[float], title: str, log_scale=False) -> None:
	spacial_image = "(cm)" in x_label and "(cm)" in y_label
	bins = []
	for label, values in [(x_label, x_values), (y_label, y_values)]:
		minimum, maximum = np.min(values), np.max(values)
		if "(%)" in label:
			bin_width = 1
		elif "(cm)" in label:
			bin_width = .03
		elif "(μm)" in label:
			bin_width = max((maximum - minimum)/80, .1)
		else:
			bin_width = (maximum - minimum)/80
		bins.append(np.arange(floor(minimum/bin_width),
		                      ceil(maximum/bin_width) + 2)*bin_width)
	counts, _, _ = np.histogram2d(x_values, y_values, bins=bins)
	vmax = np.quantile(counts, .999)
	if log_scale and vmax > 1e3:
		norm = colors.SymLogNorm(
			vmin=0, linthresh=max(30, vmax/1e3), vmax=vmax,
			linscale=1/log(10),
		)
	else:
		norm = colors.Normalize(vmin=0, vmax=vmax)
	plt.imshow(counts.T,
	           extent=(
	               np.min(x_values), np.max(x_values),
	               np.min(y_values), np.max(y_values)
	           ),
	           aspect="equal" if spacial_image else "auto",
	           norm=norm,
	           cmap=CMAP["coffee"], origin="lower")
	plt.colorbar().set_label("Counts per pixel")
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	plt.tight_layout()


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
