# analyze CPS

import os
import re
import sys
from math import log, floor, ceil, pi, sqrt
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from cr39py import cr39
from matplotlib import colors
from matplotlib.backend_bases import MouseEvent, MouseButton
from numpy.typing import NDArray
from scipy import interpolate

from cmap import CMAP

SLIT_WIDTH = .2  # (cm)

MAX_CONTRAST = 35  # (%)
MAX_ECCENTRICITY = 15  # (%)
MAX_DIAMETER = 25  # (μm)
MIN_Y = -1.2  # (cm)
MAX_Y = 0.9  # (cm)

CPS1_DISTANCE = 255  # (cm)
CPS2_DISTANCE = 255  # (cm)

def main(cps1_finger: str, cps2_finger: str, directory: str) -> None:
	for filename in os.listdir(directory):
		if filename.endswith(".cpsa"):
			# load the calibration data from disk
			calibration = load_calibration(filename, cps1_finger, cps2_finger)
			left, right = np.min(calibration.x), np.max(calibration.x)

			# load the tracks from the cpsa file
			tracks_x, tracks_y, tracks_d, tracks_c = load_tracks(directory, filename)

			# compute the spacial cuts
			data_region = (tracks_x >= left) & (tracks_x <= right) & \
			              (tracks_y >= MIN_Y) & (tracks_y <= MAX_Y)

			# ask the user about the diameter cuts
			min_diameter, max_diameter = choose_signal_region(tracks_x[data_region],
			                                                  tracks_d[data_region])
			signal = data_region & \
			         (tracks_d >= min_diameter(tracks_x)) & \
			         (tracks_d <= max_diameter(tracks_x))

			# plot the data
			plt.figure()
			histogram2d("x (cm)", tracks_x, "y (cm)", tracks_y, filename[:-4])
			plt.plot([left, left, right, right, left],
			         [MIN_Y, MAX_Y, MAX_Y, MIN_Y, MIN_Y], "k")
			plt.figure()
			histogram2d("d (μm)", tracks_d[signal], "c (%)", tracks_c[signal],
			            filename[:-4], log_scale=True)
			plt.figure()
			plt.fill_between(calibration.x,
			                 calibration.minimum_energy,
			                 calibration.maximum_energy, alpha=.5)
			plt.plot(calibration.x, calibration.nominal_energy)
			plt.xlabel("x (cm)")
			plt.ylabel("Proton energy (MeV)")

			# analyze the data
			energy, spectrum, spectrum_error = infer_spectrum(tracks_x[signal], calibration)

			# plot the results
			plt.figure()
			bar_plot("Energy (MeV)", energy, "Spectrum (MeV^-1)", spectrum, spectrum_error)
			plt.show()


def load_calibration(filename: str, cps1_finger: str, cps2_finger: str) -> "CPS":
	if "cps1" in filename.lower():
		cps = 1
	elif "cps2" in filename.lower():
		cps = 2
	elif cps2_finger.lower().startswith("n"):
		cps = 1
	elif cps1_finger.lower().startswith("n"):
		cps = 2
	else:
		raise ValueError(f"the filename doesn't make it clear whether CPS1 or CPS2 was used: "
		                 f"`{filename}`")
	if cps == 1:
		finger = cps1_finger
		slit_distance = CPS1_DISTANCE
	elif cps == 2:
		finger = cps2_finger
		slit_distance = CPS2_DISTANCE
	else:
		raise ValueError()
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
	return CPS(cps, finger, slit_distance, SLIT_WIDTH, x, energy[0, :], energy[1, :], energy[2, :])


def load_tracks(directory: str, filename: str
                ) -> tuple[NDArray[float], NDArray[float], NDArray[float], NDArray[float]]:
	file = cr39.CR39(os.path.join(directory, filename))
	file.add_cut(cr39.Cut(cmin=MAX_CONTRAST))
	file.add_cut(cr39.Cut(emin=MAX_ECCENTRICITY))
	file.add_cut(cr39.Cut(dmin=MAX_DIAMETER))
	file.apply_cuts()
	if file.ntracks == 0:
		raise ValueError("the file is now empty")
	tracks_x, tracks_y, tracks_d, tracks_c = file.trackdata_subset[:, [0, 1, 2, 3]].T
	return tracks_x, tracks_y, tracks_d, tracks_c


def choose_signal_region(x_list: NDArray[float], d_list: NDArray[float]
                         ) -> tuple[Callable[[NDArray[float]], NDArray[float]], Callable[[NDArray[float]], NDArray[float]]]:
	left, right = np.min(x_list), np.max(x_list)

	fig = plt.figure(1)
	histogram2d("x (cm)", x_list, "Diameter (μm)", d_list,
	            "click on the plot to select the minimum and maximum diameter, "
	            "then close this window.")
	lines = [plt.plot([], [], "k-")[0], plt.plot([], [], "k-")[0]]
	cursor, = plt.plot([], [], "ko")
	# if default is not None:
	# 	default_cuts, = plt.plot(default[:, 0], default[:, 1], "k-", alpha=0.3)
	# else:
	# 	default_cuts, = plt.plot([], []) TODO: load the previus one as a default

	cuts: list[list[Point]] = []

	def on_click(event: MouseEvent):
		# whenever the user clicks...
		if type(event) is MouseEvent:
			# if it's a right-click, delete a point
			if event.button == MouseButton.RIGHT:
				if len(cuts) > 0:
					if len(cuts[-1]) > 1:
						cuts[-1].pop()
					else:
						cuts.pop()
			else:
				# determine whether they are continuing a line or starting a new one
				if len(cuts) == 0 or event.xdata < cuts[-1][-1].x:
					cuts.append([])
				if len(cuts) > 2:
					raise ValueError("there should only be two lines: one above the signal and one "
					                 "below.")
				# either way, save the recent click as a new point
				cuts[-1].append(Point(event.xdata, event.ydata))
			# then update the plot
			# default_cuts.set_visible(False)
			for line in lines:
				line.set_visible(False)
			for cut, line in zip(cuts, lines):
				line.set_visible(True)
				line.set_xdata([left] + [point.x for point in cut] + [right])
				line.set_ydata([cut[0].y] + [point.y for point in cut] + [cut[-1].y])
			if len(cuts) >= 1:
				cursor.set_visible(True)
				cursor.set_xdata([cuts[-1][-1].x])
				cursor.set_ydata([cuts[-1][-1].y])
			else:
				cursor.set_visible(False)
	fig.canvas.mpl_connect('button_press_event', on_click)

	while plt.fignum_exists(1):
		plt.pause(.1)

	# once the user is done, process the results into interpolator functions
	if len(cuts) < 1:
		cuts.append([Point(0, np.max(d_list))])
	if len(cuts) < 2:
		cuts.append([Point(0, 0)])
	cuts = sorted(cuts, key=lambda line: line[0].y)
	interpolators = []
	for cut in cuts:
		cut = [Point(left, cut[0].y)] + cut + [Point(right, cut[-1].y)]
		interpolators.append(interpolate.interp1d([point.x for point in cut],
		                                          [point.y for point in cut],
		                                          bounds_error=False))
	return interpolators[0], interpolators[1]


def infer_spectrum(x_list: NDArray[float], calibration: "CPS"
                   ) -> tuple[NDArray[float], NDArray[float], NDArray[float]]:
	efficiency = (MAX_Y - MIN_Y)*calibration.slit_width/(4*pi*calibration.slit_distance**2)
	energy_bins = np.linspace(np.min(calibration.nominal_energy),
	                          np.max(calibration.nominal_energy),
	                          ceil(sqrt(x_list.size)))
	x_bins = np.interp(energy_bins, calibration.nominal_energy, calibration.x)
	counts, _ = np.histogram(x_list, x_bins)
	errors = np.sqrt(counts + 1)
	return energy_bins, counts/efficiency, errors/efficiency


def histogram2d(x_label: str, x_values: NDArray[float], y_label: str, y_values: NDArray[float],
                title: str, log_scale=False) -> None:
	spacial_image = "(cm)" in x_label and "(cm)" in y_label
	bins = []
	for label, values in [(x_label, x_values), (y_label, y_values)]:
		minimum, maximum = np.min(values), np.max(values)
		if "(%)" in label:
			bin_width, bin_number, quantized = 1, None, True
		elif "(cm)" in label:
			bin_width, bin_number, quantized = .03, None, False
		elif "(μm)" in label:
			bin_width, bin_number, quantized = None, min(80, floor((maximum - minimum)/.1)), False
		else:
			bin_width, bin_number, quantized = None, 80, False
		if quantized:
			bins.append(np.arange(floor(minimum/bin_width),
			                      floor(maximum/bin_width) + 2)*bin_width)
		else:
			if bin_number is None:
				bin_number = round((maximum - minimum)/bin_width)
			bins.append(np.linspace(minimum, maximum, bin_number + 1))
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


def bar_plot(x_label: str, bar_edges: NDArray[float],
             y_label: str, bar_heights: NDArray[float], bar_errors: NDArray[float]) -> None:
	x = np.repeat(bar_edges, 2)[1:-1]
	y = np.repeat(bar_heights, 2)
	plt.plot(x, y, "k-", linewidth=1)
	bar_centers = (bar_edges[:-1] + bar_edges[1:])/2
	plt.errorbar(x=bar_centers, y=bar_heights, yerr=bar_errors, fmt="k-", linewidth=1)
	plt.xlabel(x_label)
	plt.ylabel(y_label)


class Point:
	def __init__(self, x: float, y: float):
		self.x = x
		self.y = y


class CPS:
	def __init__(self, cps: int, finger: str, slit_distance: float, slit_width: float,
	             x: NDArray[float], minimum_energy: NDArray[float],
	             nominal_energy: NDArray[float], maximum_energy: NDArray[float]):
		self.cps = cps
		self.finger = finger
		self.slit_distance = slit_distance
		self.slit_width = slit_width
		self.x = x
		self.minimum_energy = minimum_energy
		self.nominal_energy = nominal_energy
		self.maximum_energy = maximum_energy


if __name__ == "__main__":
	if len(sys.argv) == 4:
		cps1_finger, cps2_finger, directory = sys.argv[1:]
		with open("arguments.txt", "w") as file:
			file.write(f"cps1-finger={cps1_finger}\n"
			           f"cps2-finger={cps2_finger}\n"
			           f"directory={directory}\n")
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
							raise ValueError(f"The `arguments.txt` file contains an unrecognized "
							                 f"key: {key}")
		except IOError:
			raise ValueError("You must run this script with three command line arguments: "
			                 "`python analyze_cps.py cps1-finger cps2-finger directory`; "
			                 "or specify the arguments by creating a file called `arguments.txt` "
			                 "formatted like\n"
			                 "  “cps1-finger=a1\n   cps2-finger=b2w\n   directory=example/path/”")
		if cps1_finger is None or cps2_finger is None or directory is None:
			raise ValueError("The `arguments.txt` file was missing some of the three required "
			                 "arguments.")
	else:
		raise ValueError("You must run this script with exactly three command line arguments: "
		                 "`python analyze_cps.py cps1-finger cps2-finger directory`, "
		                 "or alternatively specify the arguments by creating a file called "
		                 "`arguments.txt` formatted like \n"
		                 "  “cps1-finger=a1\n   cps2-finger=b2w\n   directory=example/path/”")

	main(cps1_finger, cps2_finger, directory)
