# analyze CPS

import os
import re
import sys
from math import log, floor, ceil, pi, sqrt
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cr39py import cr39
from matplotlib import colors
from matplotlib.backend_bases import MouseEvent, MouseButton
from numpy.typing import NDArray
from pandas import DataFrame
from scipy import interpolate
from xarray import DataArray

from analyze_spectra import plot_bars, Spectrum, FIGURE_SIZE
from cmap import CMAP

SLIT_WIDTH = 2e-1  # (cm)

DATA_REGION = (-1.2, 0.9)  # y_min, y_max (cm)

MAX_CONTRAST = 20  # (%)
MAX_ECCENTRICITY = 15  # (%)
MAX_DIAMETER = 25  # (μm)
BIN_SIZE = .05  # (MeV)

CPS1_DISTANCE = 255  # (cm)
CPS2_DISTANCE = 100  # (cm)

X = "x (cm)"
Y = "y (cm)"
D = "Track diameter (μm)"
C = "Track contrast (%)"
SPACIAL_DIMS = {X, Y}


def main(cps1_finger: str, cps2_finger: str, particle: str, directory: str) -> None:
	for filename in os.listdir(directory):
		if filename.endswith(".cpsa"):
			# parse the particle type
			particle_name, particle_mass = parse_particle(particle)

			# load the calibration data from disk
			calibration = load_calibration(filename, cps1_finger, cps2_finger, particle_mass)
			left, right = np.min(calibration.x), np.max(calibration.x)

			# load the tracks from the cpsa file
			tracks = load_tracks(os.path.join(directory, filename))
			data = in_rectangle(tracks, left, right, *DATA_REGION)

			# ask the user about the background region
			background_region = choose_background_region(tracks, left, right, *DATA_REGION)
			background = calculate_background(tracks, *background_region)

			# ask the user about the diameter cuts
			min_diameter_cut, max_diameter_cut = choose_diameter_cuts(tracks[data], background)
			signal = data & apply_diagonal_cuts(tracks, min_diameter_cut, max_diameter_cut)

			# plot the cleaned-up data
			plt.figure(figsize=FIGURE_SIZE)
			plot_2d_histogram(tracks, X, Y, filename[:-5])
			plot_rectangle(*background_region, label="Background")
			plot_rectangle(left, right, *DATA_REGION, label="Data")
			plt.tight_layout()
			plt.figure(figsize=FIGURE_SIZE)
			plot_2d_histogram(tracks[signal], D, C, filename[:-5], background, log_scale=True)
			plt.tight_layout()
			plt.figure(figsize=FIGURE_SIZE)
			plt.fill_between(calibration.x,
			                 calibration.minimum_energy,
			                 calibration.maximum_energy, alpha=.5)
			plt.plot(calibration.x, calibration.nominal_energy)
			plt.xlabel(X)
			plt.ylabel(f"{particle_name} energy (MeV)")
			plt.tight_layout()

			# analyze and save the data
			spectrum = infer_spectrum(
				tracks[signal], calibration, background, min_diameter_cut, max_diameter_cut)
			spectrum = downsample(spectrum)
			save_spectrum(spectrum, particle_name, directory, filename[:-5])

			# plot the results
			plt.figure(figsize=FIGURE_SIZE)
			plot_bars(spectrum, f"{particle_name} energy (MeV)", "Spectrum (MeV^-1)")
			plt.tight_layout()

			plt.show()


def load_calibration(filename: str, cps1_finger: str, cps2_finger: str, particle_mass: float) -> "CPS":
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

	energy /= particle_mass
	return CPS(cps, finger, slit_distance, SLIT_WIDTH, x, energy[0, :], energy[1, :], energy[2, :])


def load_tracks(filepath: str) -> DataFrame:
	file = cr39.CR39(filepath)
	file.add_cut(cr39.Cut(cmin=MAX_CONTRAST))
	file.add_cut(cr39.Cut(emin=MAX_ECCENTRICITY))
	file.add_cut(cr39.Cut(dmin=MAX_DIAMETER))
	file.apply_cuts()
	if file.ntracks == 0:
		raise ValueError("the file is empty")
	headers = [X, Y, D, C]
	columns = {header: file.trackdata_subset[:, i] for i, header in enumerate(headers)}
	return DataFrame(columns)


def parse_particle(code: str) -> tuple[str, float]:
	if code.lower().startswith("p"):
		mass = 1
	elif code.lower().startswith("d"):
		mass = 2
	elif code.lower().startswith("t"):
		mass = 3
	elif code.lower().startswith("a"):
		mass = 1/2
	else:
		try:
			mass = float(code)
		except ValueError:
			raise ValueError(f"Unrecognized charged particle: '{code}'")

	if round(mass, 1) == 0.5:
		name = "Alpha"
	elif round(mass, 1) == 1.0:
		name = "Proton"
	elif round(mass, 1) == 2.0:
		name = "Deuteron"
	elif round(mass, 1) == 3.0:
		name = "Triton"
	else:
		name = "Particle"
	return name, mass


def calculate_background(data: DataFrame, left: float, right: float, bottom: float, top: float) -> DataArray:
	data = data[in_rectangle(data, left, right, bottom, top)]
	d_bin_edges = get_bin_edges(D, data[D])
	c_bin_edges = get_bin_edges(C, data[C])
	counts, _, _ = np.histogram2d(data[D], data[C], bins=(d_bin_edges, c_bin_edges))
	counts = DataArray(counts, dims=(D, C), coords={D: d_bin_edges[0:-1], C: c_bin_edges[0:-1]})
	area = (right - left) * (top - bottom)
	return counts/area


def in_rectangle(data: DataFrame,
                 left: float, right: float, bottom: float, top: float) -> NDArray[bool]:
	return (data[X] >= left) & (data[X] <= right) & \
	       (data[Y] >= bottom) & (data[Y] <= top)


def choose_background_region(tracks: DataFrame, data_left: float, data_right: float,
                             data_bottom: float, data_top: float) -> tuple[float, float, float, float]:
	fig = plt.figure("selection", figsize=FIGURE_SIZE)
	plot_2d_histogram(tracks, X, Y, "click to set the corners of the background region, then close this plot")
	plot_rectangle(data_left, data_right, data_bottom, data_top, label="Data region")
	plt.tight_layout()
	points, = plt.plot([], [], "ko")
	rectangle, = plt.plot([], [], "k-")

	vertices: list[Point] = []

	def on_click(event: MouseEvent):
		# whenever the user clicks...
		if type(event) is MouseEvent and event.xdata is not None:
			# if it's a right-click, delete a point
			if event.button == MouseButton.RIGHT:
				if len(vertices) > 0:
					vertices.pop()
			# otherwise, save a new point
			elif len(vertices) < 2:
				vertices.append(Point(event.xdata, max(event.ydata, DATA_REGION[1] + .1)))
			# then update the plot
			points.set_xdata([vertex.x for vertex in vertices])
			points.set_ydata([vertex.y for vertex in vertices])
			if len(vertices) >= 2:
				rectangle.set_visible(True)
				a, b = vertices[:2]
				rectangle.set_xdata([a.x, b.x, b.x, a.x, a.x])
				rectangle.set_ydata([a.y, a.y, b.y, b.y, a.y])
			else:
				rectangle.set_visible(False)
	fig.canvas.mpl_connect('button_press_event', on_click)

	while plt.fignum_exists("selection"):
		plt.pause(.1)
	if len(vertices) != 2:
		raise ValueError("you didn't specify both corners of the rectangle.")

	# once the user is done, arrange the results
	a, b = vertices
	left, right = min(a.x, b.x), max(a.x, b.x)
	bottom, top = min(a.y, b.y), max(a.y, b.y)
	return left, right, bottom, top


def choose_diameter_cuts(tracks: DataFrame, background: DataArray,
                         ) -> tuple[list["Point"], list["Point"]]:
	left, right = tracks[X].min(), tracks[X].max()

	fig = plt.figure("selection", figsize=FIGURE_SIZE)
	plot_2d_histogram(tracks, X, D,
	                  "click on the plot to select the minimum and maximum diameter, "
	                  "then close this window.", background)
	plt.tight_layout()
	lines = [plt.plot([], [], "k-")[0], plt.plot([], [], "k-")[0]]
	cursor, = plt.plot([], [], "ko")
	# if default is not None:
	# 	default_cuts, = plt.plot(default[:, 0], default[:, 1], "k-", alpha=0.3)
	# else:
	# 	default_cuts, = plt.plot([], []) TODO: load the previus one as a default

	cuts: list[list[Point]] = []

	def on_click(event: MouseEvent):
		# whenever the user clicks...
		if type(event) is MouseEvent and event.xdata is not None:
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
					if len(cuts) < 2:
						cuts.append([])
					else:
						event.xdata = cuts[-1][-1].x
				# either way, save the recent click as a new point
				cuts[-1].append(Point(event.xdata, event.ydata))
			# then update the plot
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

	while plt.fignum_exists("selection"):
		plt.pause(.1)

	# once the user is done, process the results into interpolator functions
	if len(cuts) < 1:
		cuts.append([Point(0, tracks[D].max())])
	if len(cuts) < 2:
		cuts.append([Point(0, 0)])
	cuts = sorted(cuts, key=lambda line: line[0].y)
	for cut in cuts:
		cut.insert(0, Point(left, cut[0].y))
		cut.append(Point(right, cut[-1].y))
	return cuts[0], cuts[1]


def apply_diagonal_cuts(data: DataFrame, minimum_diameter: list["Point"],
                        maximum_diameter: list["Point"]) -> NDArray[bool]:
	minimum_diameter_at = interpolate.interp1d([p.x for p in minimum_diameter],
	                                           [p.y for p in minimum_diameter], bounds_error=False)
	maximum_diameter_at = interpolate.interp1d([p.x for p in maximum_diameter],
	                                           [p.y for p in maximum_diameter], bounds_error=False)
	return (data[D] >= minimum_diameter_at(data[X])) & (data[D] <= maximum_diameter_at(data[X]))


def get_bin_edges(label: str, values: NDArray[float]) -> DataArray:
	minimum, maximum = np.min(values), np.max(values)
	if "(%)" in label:
		bin_width, num_bins, quantized = 1, None, True
	elif "(cm)" in label:
		bin_width, num_bins, quantized = .03, None, False
	elif "(μm)" in label:
		bin_width, num_bins, quantized = None, min(80, floor((maximum - minimum)/.1)), False
	else:
		bin_width, num_bins, quantized = None, 80, False
	if quantized:
		return np.arange(-0.5, maximum/bin_width + 1)*bin_width
	else:
		if num_bins is None:
			num_bins = round((maximum - minimum)/bin_width)
		return DataArray(np.linspace(minimum, maximum, num_bins + 1), dims=(label,))


def infer_spectrum(data: DataFrame, calibration: "CPS", background: DataArray,
                   min_diameter_cut: list["Point"], max_diameter_cut: list["Point"],
                   ) -> "Spectrum":
	# calculate the scalar prefactor
	slit_height = DATA_REGION[1] - DATA_REGION[0]
	efficiency = slit_height*calibration.slit_width/(4*pi*calibration.slit_distance**2)

	# compute the x bins by converting from energy bins
	energy_bin_edges = np.linspace(np.min(calibration.nominal_energy),
	                               np.max(calibration.nominal_energy),
	                               ceil(np.ptp(calibration.nominal_energy)/BIN_SIZE))
	energy_bin_width = energy_bin_edges[1] - energy_bin_edges[0]
	x_bin_edges = DataArray(
		np.interp(energy_bin_edges, calibration.nominal_energy, calibration.x), dims=(X,))

	# do the histogramming using the x bins
	counts = DataArray(np.histogram(data[X], x_bin_edges)[0], dims=(X,))
	errors = np.sqrt(counts + 1)  # TODO: this doesn't account for uncertainty in the background

	# arrange the background array's dimensions as needed and do the background subtraction
	for dim in background.dims:
		if dim != D:
			background = background.sum(dim=dim)
	for dim in SPACIAL_DIMS:
		if dim != X:
			background = background*(data[dim].max() - data[dim].min())  # TODO: it would be better to store the range information from when we made the cuts
	background = background*(x_bin_edges[1:] - x_bin_edges[0:-1])
	d_bins = background.coords[D]
	#  do a little numerical integral to see how much diameter from each bin falls within the d cuts
	binned_background = np.zeros(counts.shape)
	offsets = np.arange(0.5, 6)/6
	for dx in offsets:
		for dd in offsets:
			d = d_bins + dd*(d_bins[1] - d_bins[0])
			x = x_bin_edges[0:-1] + dx*(x_bin_edges[1:] - x_bin_edges[0:-1])
			d_min = DataArray(
				np.interp(x, [p.x for p in min_diameter_cut], [p.y for p in min_diameter_cut]), dims=(X,))
			d_max = DataArray(
				np.interp(x, [p.x for p in max_diameter_cut], [p.y for p in max_diameter_cut]), dims=(X,))
			signal = (d >= d_min) & (d <= d_max)
			binned_background += xr.where(signal, background, 0).sum(dim=D)/offsets.size**2
	counts = counts - binned_background

	spectral_density = counts/energy_bin_width/efficiency
	spectral_error = errors/energy_bin_width/efficiency

	return Spectrum(energy_bin_edges, spectral_density.to_numpy(), spectral_error.to_numpy())


def downsample(spectrum: "Spectrum") -> "Spectrum":
	typical_value = np.quantile(spectrum.values, .90)
	typical_error = np.quantile(spectrum.errors, .90)
	factor = max(1, round(sqrt(typical_error/typical_value/.05)))
	indices = np.reshape(np.arange(floor(spectrum.values.size/factor)*factor), (-1, factor))
	return Spectrum(spectrum.energy_bin_edges[0::factor],
	                spectrum.values[indices].mean(axis=1),
	                (spectrum.errors[indices]**-2).sum(axis=1)**(-1/2))


def plot_rectangle(left: float, right: float, bottom: float, top: float, *,
                   label: Optional[str] = None) -> None:
	plt.plot([left, right, right, left, left],
	         [bottom, bottom, top, top, bottom], "k")
	if label is not None:
		plt.text((left + right)/2, (bottom + top)/2, label)


def plot_2d_histogram(data: DataFrame, x_label: str, y_label: str, title: str,
                      background: Optional[DataArray] = None, log_scale=False) -> None:
	# set up the binning
	spacial_image = x_label in SPACIAL_DIMS and y_label in SPACIAL_DIMS
	x_bin_edges = get_bin_edges(x_label, data[x_label])
	y_bin_edges = get_bin_edges(y_label, data[y_label])

	# compute the histogram
	counts = DataArray(
		np.histogram2d(data[x_label], data[y_label], bins=(x_bin_edges, y_bin_edges))[0],
		dims=(x_label, y_label))

	# subtract the background...
	if background is not None:
		for dim in background.dims:
			if dim != x_label and dim != y_label:
				background = background.sum(dim=dim)
		for dim in SPACIAL_DIMS:
			if dim == x_label:
				background = background*(x_bin_edges[1] - x_bin_edges[0])
			elif dim == y_label:
				background = background*(y_bin_edges[1] - y_bin_edges[0])
			else:
				background = background*(data[dim].max() - data[dim].min())
		counts -= background

	# set up the limits
	vmax = np.quantile(counts, .999)
	if log_scale and vmax > 1e3:
		norm = colors.SymLogNorm(
			vmin=0, linthresh=max(30, vmax/1e3), vmax=vmax,
			linscale=1/log(10),
		)
	else:
		norm = colors.Normalize(vmin=0, vmax=vmax)

	# make the plot
	plt.imshow(counts.T,
	           extent=(
		           data[x_label].min(), data[x_label].max(),
		           data[y_label].min(), data[y_label].max(),
	           ),
	           aspect="equal" if spacial_image else "auto",
	           norm=norm,
	           cmap=CMAP["coffee"], origin="lower")
	plt.colorbar().set_label("Counts per pixel")
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	plt.tight_layout()


def save_spectrum(spectrum: "Spectrum",
                  particle: str, directory, filename: str) -> None:
	energies = (spectrum.energy_bin_edges[0:-1] + spectrum.energy_bin_edges[1:])/2
	dataframe = DataFrame({f"{particle} energy (MeV)": energies,
	                       "Spectrum (MeV^-1)": spectrum.values,
	                       "Spectrum error (MeV^-1)": spectrum.errors})
	dataframe.to_csv(os.path.join(directory, filename + "_spectrum.csv"), index=False)


class Point:
	def __init__(self, x: float, d: float):
		self.x = x
		self.y = d

	def __str__(self) -> str:
		return f"Point({self.x}, {self.y})"


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

	def __str__(self) -> str:
		return f"CPS({self.cps}, {self.finger}, {self.slit_distance}, {self.slit_width}, ...)"


if __name__ == "__main__":
	if len(sys.argv) == 5:
		cps1_finger, cps2_finger, particle, directory = sys.argv[1:]
		with open("arguments.txt", "w") as file:
			file.write(f"cps1-finger={cps1_finger}\n"
			           f"cps2-finger={cps2_finger}\n"
			           f"particle={particle}\n"
			           f"directory={directory}\n")
	elif len(sys.argv) == 1:
		cps1_finger, cps2_finger, particle, directory = None, None, None, None
		try:
			with open("arguments.txt", "r") as file:
				for line in file:
					if "=" in line:
						key = line[:line.index("=")].lower()
						value = line[line.index("=") + 1:].lower().strip()
						if "cps1" in key:
							cps1_finger = value
						elif "cps2" in key:
							cps2_finger = value
						elif "particle" in key or "ion" in key or "a/z^2" in key:
							particle = value
						elif "directory" in key or "path" in key:
							directory = value
						else:
							raise ValueError(f"The `arguments.txt` file contains an unrecognized "
							                 f"key: {key}")
		except IOError:
			raise ValueError("You must run this script with four command line arguments: "
			                 "`python analyze_cps.py cps1-finger cps2-finger directory`; "
			                 "or specify the arguments by creating a file called `arguments.txt` "
			                 "formatted like\n"
			                 "  “cps1-finger=a1\n   cps2-finger=b2w\n"
			                 "   particle=d\n   directory=example/path/”")
		if cps1_finger is None or cps2_finger is None or particle is None or directory is None:
			raise ValueError("The `arguments.txt` file was missing some of the four required "
			                 "arguments: cps1, cps2, key, and directory.")
	else:
		raise ValueError("You must run this script with exactly four command line arguments: "
		                 "`python analyze_cps.py cps1-finger cps2-finger directory`, "
		                 "or alternatively specify the arguments by creating a file called "
		                 "`arguments.txt` formatted like \n"
		                 "  “cps1-finger=a1\n   cps2-finger=b2w\n"
		                 "   particle=d\n   directory=example/path/”")

	main(cps1_finger, cps2_finger, particle, directory)
