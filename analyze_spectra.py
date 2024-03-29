import argparse
import os
import re
from math import sqrt, floor, pi, log, log10, nan, inf
from typing import Union, Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backend_bases import MouseEvent, MouseButton
from numpy.typing import NDArray
from pandas.errors import ParserError
from scipy import optimize

FIGURE_SIZE = (7.5, 4)
HIGHLIGHT_COLOR = "#C00000"
CATEGORICAL_COLORS = ["#406836", "#BA5662", "#1D4881"]


def main(directory: str) -> None:
	# load all the spectra and organize them
	spectra: dict[str, dict[int, list[tuple[Spectrum, str, str, str]]]] = {}
	for subdirectory, _, filenames in os.walk(directory):
		for filename in filenames:
			if filename.endswith(".csv"):
				cps_descriptor = re.search(r"cps(\d+)_?", filename, re.IGNORECASE)
				if cps_descriptor is not None:
					cps = int(cps_descriptor.group(1))
				else:
					cps = 0
				shot_descriptor = re.search(r"O?\d{5,6}", filename, re.IGNORECASE)
				if shot_descriptor is not None:
					shot = shot_descriptor.group(0)
				else:
					shot = filename[:-4]
				if shot not in spectra:
					spectra[shot] = {}
				if cps not in spectra[shot]:
					spectra[shot][cps] = []
				try:
					spectrum, energy_label, spectrum_label = load_spectrum(
						os.path.join(subdirectory, filename))
				except (ParserError, InvalidFileError):
					print(f"I'm skipping {filename} because I can't read it.")
					continue
				spectra[shot][cps].append((spectrum, energy_label, spectrum_label, filename[:-4]))
	if len(spectra) == 0:
		print(f"No CPS spectra were found in `{directory}`.")

	# then go thru them by shot
	for shot in sorted(spectra.keys()):
		for spectra_on_each_cps in spectra[shot].values():
			for spectrum, energy_label, spectrum_label, filename in spectra_on_each_cps:
				print(filename)

				# ask the user for the energy limits
				left, right = choose_limits(spectrum, energy_label, spectrum_label)

				# get the important numbers
				raw = count_raw_statistics(spectrum, left, right)
				try:
					gaussian = fit_gaussian(spectrum, left, right)
				except RuntimeError as e:
					print(e)
					continue

				# plot the results
				plt.figure(figsize=FIGURE_SIZE)
				plt.locator_params(steps=[1, 2, 5, 10])
				plot_bars(spectrum, energy_label, spectrum_label)
				plt.plot(spectrum.energy_bin_edges,
				         gaussian_function(spectrum.energy_bin_edges, gaussian),
				         HIGHLIGHT_COLOR, zorder=2)
				plt.title(filename)
				annotate_plot(f"Raw yield = {raw.total:.2e}\n"
				              f"Raw mean = {raw.mean:.2f} MeV\n"
				              f"Fit yield = {gaussian.total:.2e}\n"
				              f"Fit mean = {gaussian.mean:.2f} MeV\n"
				              f"Fit width = {gaussian.sigma*2*sqrt(2*log(2)):.2f} MeV")
				plt.tight_layout()

				# save and display the figure
				plt.savefig(os.path.join(directory, filename + ".png"),
				            dpi=300, transparent=True)
				plt.show()

		# add in an overlaid plot if both CPS were used
		if len(spectra[shot]) > 1:
			# do one for each set of corresponding fingers
			num_finger_pairs = min(len(spectra[shot][cps]) for cps in spectra[shot])
			for finger_index in range(num_finger_pairs):
				plt.figure(figsize=FIGURE_SIZE)
				plt.locator_params(steps=[1, 2, 5, 10])
				energy_minima, energy_maxima = [], []
				for cps, spectra_on_this_cps in spectra[shot].items():
					spectrum, energy_label, spectrum_label, _ = spectra_on_this_cps[finger_index]
					plot_bars(spectrum, energy_label, spectrum_label,
					          color=CATEGORICAL_COLORS[cps], label=f"CPS{cps}")
					energy_minima.append(spectrum.energy_bin_edges[0])
					energy_maxima.append(spectrum.energy_bin_edges[-1])
				plt.legend()
				plt.xlim(max(energy_minima), min(energy_maxima))
				plt.tight_layout()
				plt.savefig(os.path.join(directory, f"{shot}_{finger_index}_spectra.png"), transparent=True)

		plt.show()


def load_spectrum(filepath: str) -> tuple["Spectrum", str, str]:
	""" load a spectrum previusly saved as a CSV file """
	data = pd.read_csv(filepath)
	energy_bins, values, errors = None, None, None
	energy_label, spectrum_label = None, None
	for column in data.columns:
		if "energy" in column.lower():
			energy_bins = data[column].to_numpy()
			energy_label = column
		elif "error" in column.lower():
			errors = data[column].to_numpy()
		elif "spectrum" in column.lower():
			values = data[column].to_numpy()
			spectrum_label = column

	if energy_bins is None or values is None or errors is None:
		raise InvalidFileError("I found a file that didn't have the right columns")

	bin_width = energy_bins[1] - energy_bins[0]
	energy_bin_edges = np.concatenate([
		energy_bins - bin_width/2,
		[energy_bins[-1] + bin_width/2]])

	return Spectrum(energy_bin_edges, values, errors), energy_label, spectrum_label


def choose_limits(spectrum: "Spectrum", x_label: str, y_label: str) -> tuple[float, float]:
	""" prompt the user to click on a plot to choose lower and upper limits """
	fig = plt.figure("selection", figsize=FIGURE_SIZE)
	plt.locator_params(steps=[1, 2, 5, 10])
	plot_bars(spectrum, x_label, y_label)
	plt.title("click to select the lower and upper bounds of the peak, then close this plot")
	plt.tight_layout()
	lines = [plt.plot([], [], "k--")[0], plt.plot([], [], "k--")[0]]
	curve, = plt.plot(spectrum.energy_bin_edges, np.zeros_like(spectrum.energy_bin_edges),
	                  HIGHLIGHT_COLOR, zorder=2)
	curve.set_visible(False)

	minimum = spectrum.energy_bin_edges[0]
	maximum = spectrum.energy_bin_edges[-1]
	limits: list[float] = []

	def update_plot():
		nonlocal minimum, maximum
		if len(limits) < 2:
			minimum = spectrum.energy_bin_edges[0]
			maximum = spectrum.energy_bin_edges[-1]
		else:
			minimum, maximum = min(limits), max(limits)
		try:
			gaussian = fit_gaussian(spectrum, minimum, maximum)
		except RuntimeError:
			curve.set_visible(False)
		else:
			curve.set_ydata(gaussian_function(spectrum.energy_bin_edges, gaussian))
			curve.set_visible(True)
	update_plot()

	def on_click(event: MouseEvent):
		# whenever the user clicks...
		if type(event) is MouseEvent and event.xdata is not None:
			# if it's a right-click, delete a point
			if event.button == MouseButton.RIGHT:
				if len(limits) > 0:
					limits.pop()
			# otherwise, save a new point
			elif len(limits) < 2:
				limits.append(event.xdata)
			# then update the plot
			for i in range(len(lines)):
				if i < len(limits):
					lines[i].set_xdata([limits[i], limits[i]])
					lines[i].set_ydata([np.min(spectrum.values), np.max(spectrum.values)])
					lines[i].set_visible(True)
				else:
					lines[i].set_visible(False)
			update_plot()
	fig.canvas.mpl_connect('button_press_event', on_click)

	while plt.fignum_exists("selection"):
		plt.pause(.1)

	# once the user is done, arrange the results
	return minimum, maximum


def count_raw_statistics(spectrum: "Spectrum", left: float, right: float) -> "Distribution":
	# first convert this to a series of bin centers and bin counts
	energies = (spectrum.energy_bin_edges[0:-1] + spectrum.energy_bin_edges[1:])/2
	counts = np.where((energies >= left) & (energies <= right),
	                  spectrum.values*spectrum.energy_bin_widths,
	                  0)
	errors = spectrum.errors*spectrum.energy_bin_widths
	# then do the math
	total = Quantity(
		np.sum(counts),  # type: ignore
		sqrt(np.sum(errors**2)))
	mean = Quantity(
		np.sum(energies*counts)/total.value,
		sqrt(np.sum((errors*(energies/total.value - np.sum(energies*counts)/total.value**2))**2)))
	sigma = Quantity(
		sqrt(np.sum((energies - mean.value)**2)),
		nan)  # eh, no one’s going to look at the sigma error
	# tie it up in an object and return
	return Distribution(total, mean, sigma)


def fit_gaussian(spectrum: "Spectrum", left: float, right: float) -> "Distribution":
	""" fit a gaussian to a given spectrum within some energy bounds """
	energy_bin_centers = (spectrum.energy_bin_edges[0:-1] + spectrum.energy_bin_edges[1:])/2
	in_limits = (energy_bin_centers >= left) & (energy_bin_centers <= right)
	raw_total = np.sum(spectrum.values*spectrum.energy_bin_widths, where=in_limits)
	try:
		popt, pcov = optimize.curve_fit(f=gaussian_function,
		                                xdata=energy_bin_centers[in_limits],
		                                ydata=spectrum.values[in_limits],
		                                sigma=spectrum.errors[in_limits],
		                                p0=[raw_total, (left + right)/2, (right - left)/2])
	except (RuntimeError, ValueError, TypeError):
		raise RuntimeError("Could not find optimal Gaussian parameters!")
	values = []
	for i in range(len(popt)):
		values.append(Quantity(popt[i], sqrt(pcov[i, i])))
	return Distribution(values[0], values[1], values[2])


def gaussian_function(x: NDArray[float],
                      *args: Union["Distribution", float]) -> NDArray[float]:
	""" return the value of a gaussian curve at the given x, specifying the gaussian’s parameters """
	if len(args) == 1:
		N, μ, σ = args[0].total.value, args[0].mean.value, args[0].sigma.value
	elif len(args) == 3:
		N, μ, σ = args
	else:
		raise TypeError("gaussian_function() must be passed a Gaussian object or three floats.")
	return N/sqrt(2*pi*σ**2)*np.exp(-(x - μ)**2/(2*σ**2))


def plot_bars(spectrum: "Spectrum", x_label: str, y_label: str,
              color: str = "k", label: Optional[str] = None) -> None:
	""" plot and label a spectrum in that blocky style that makes it look like a histogram, with error bars """
	x = np.repeat(spectrum.energy_bin_edges, 2)[1:-1]
	y = np.repeat(spectrum.values, 2)
	plt.plot(x, y, color, linewidth=0.7, zorder=3, label=label)
	energy_bin_centers = (spectrum.energy_bin_edges[:-1] + spectrum.energy_bin_edges[1:])/2
	plt.errorbar(x=energy_bin_centers, y=spectrum.values, yerr=spectrum.errors,
	             ecolor=color, elinewidth=0.7, fmt="none")
	plt.xlim(spectrum.energy_bin_edges[0], spectrum.energy_bin_edges[-1])
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid("on")


def annotate_plot(text: str) -> None:
	""" put some text in the upper right-hand corner of the plot """
	text = plt.text(.99, .98, text, zorder=4,
	                ha='right', va='top', transform=plt.gca().transAxes)
	text.set_bbox(dict(facecolor='w', alpha=0.5, edgecolor="none"))


class InvalidFileError(Exception):
	""" an error thrown when the specified file is a valid CSV, but not a valid spectrum """
	pass


class Spectrum:
	def __init__(self, energy_bin_edges: NDArray[float], values: NDArray[float], errors: NDArray[float]):
		""" an object encapsulating the energy bins, counts, and error bars of a spectrum """
		assert energy_bin_edges.ndim == 1 and values.ndim == 1 and errors.ndim == 1
		assert energy_bin_edges.size - 1 == values.size and values.size == errors.size
		self.energy_bin_edges = energy_bin_edges
		self.values = values
		self.errors = errors
		self.energy_bin_widths = energy_bin_edges[1:] - energy_bin_edges[0:-1]


class Quantity:
	def __init__(self, value: float, error: float):
		""" an inferred number and its error bar """
		self.value = value
		self.error = error

	def __mul__(self, other: float) -> "Quantity":
		return Quantity(self.value*other, self.error*other)

	def __format__(self, format_spec: str) -> str:
		# for "e", manage the exponent manually to make the value and error match
		if "e" in format_spec:
			exponent = floor(log10(abs(self.value)))
			new_format_spec = format_spec.replace("e", "f")
			return f"{format(self.value/10**exponent, new_format_spec)}e{exponent:+03d} ± " \
			       f"{format(self.error/10**exponent, new_format_spec)}e{exponent:+03d}"
		# otherwise just delegate tothe built-in format for each of the value and error
		else:
			return f"{format(self.value, format_spec)} ± " \
			       f"{format(self.error, format_spec)}"


class Distribution:
	def __init__(self, total: Quantity, mean: Quantity, sigma: Quantity):
		""" the three values that define a gaussian curve """
		self.total = total
		self.mean = mean
		self.sigma = sigma


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		prog="python analyze_spectra.py",
		description="Convert analyzed spectra to useful plots and numbers.")
	parser.add_argument("--directory", type=str, default="./",
	                    help="Absolute or relative path to the folder containing the spectrum CSV files (not necessary if the scan files are located somewhere in the current working directory)")
	args = parser.parse_args()

	main(args.directory)
