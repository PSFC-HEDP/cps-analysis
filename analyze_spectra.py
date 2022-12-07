import os
import re
import sys
from math import sqrt, floor, pi, log, log10
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


def main(directory: str) -> None:
	# load all the spectra and organize them
	spectra: dict[str, dict[int, tuple[Spectrum, str, str, str]]] = {}
	for filename in os.listdir(directory):
		if filename.endswith(".csv"):
			cps_descriptor = re.search(r"cps(\d+)_?", filename.lower())
			if cps_descriptor is not None:
				cps = int(cps_descriptor.group(1))
				shot = filename.lower().replace(cps_descriptor.group(0), "")
			else:
				cps = 0
				shot = filename
			if shot not in spectra:
				spectra[shot] = {}
			try:
				spectrum, energy_label, spectrum_label = load_spectrum(
					os.path.join(directory, filename))
			except (ParserError, InvalidFileError):
				print(f"I'm skipping {filename} because I can't read it.")
				continue
			spectra[shot][cps] = downsample(spectrum), energy_label, spectrum_label, filename[:-4]

	# then go thru them by shot
	for shot in sorted(spectra.keys()):
		for spectrum, energy_label, spectrum_label, filename in spectra[shot].values():
			# ask the user for the energy limits
			left, right = choose_limits(spectrum, energy_label, spectrum_label)

			# get the important numbers
			total_yield = Quantity(np.sum(spectrum.values*spectrum.energy_bin_widths),
			                       sqrt(np.sum(spectrum.errors**2*spectrum.energy_bin_widths)))
			gaussian = fit_gaussian(spectrum, left, right)

			# plot the results
			plt.figure(figsize=FIGURE_SIZE)
			plot_bars(spectrum, energy_label, spectrum_label)
			plt.plot(spectrum.energy_bin_edges,
			         gaussian_function(spectrum.energy_bin_edges, gaussian),
			         HIGHLIGHT_COLOR, zorder=1)
			annotate_plot(f"Total yield = {total_yield:.2e}\n"
			              f"Peak yield = {gaussian.total:.2e}\n"
			              f"Peak energy = {gaussian.mean:.2f} MeV\n"
			              f"Peak width = {gaussian.sigma*2*sqrt(2*log(2)):.2f} MeV")
			plt.tight_layout()

			# save the data, and also the most recent figure
			plt.savefig(os.path.join(directory, filename + ".png"),
			            dpi=300, transparent=True)

		# add in an overlaid plot if there are a few of these
		if len(spectra[shot]) > 1:
			plt.figure(figsize=FIGURE_SIZE)
			energy_minima, energy_maxima = [], []
			for cps, (spectrum, energy_label, spectrum_label, _) in spectra[shot].items():
				plot_bars(spectrum, energy_label, spectrum_label,
				          color=f"C{cps}", label=f"CPS{cps}")
				energy_minima.append(spectrum.energy_bin_edges[0])
				energy_maxima.append(spectrum.energy_bin_edges[-1])
			plt.xlim(min(energy_minima), max(energy_maxima))
			plt.tight_layout()
			plt.savefig(os.path.join(directory, shot + "_spectra.png"))

		plt.show()


def load_spectrum(filepath: str) -> tuple["Spectrum", str, str]:
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


def downsample(spectrum: "Spectrum") -> "Spectrum":
	typical_value = np.quantile(spectrum.values, .90)
	typical_error = np.quantile(spectrum.errors, .90)
	factor = max(1, round(sqrt(typical_error/typical_value/.05)))
	indices = np.reshape(np.arange(floor(spectrum.values.size/factor)*factor), (-1, factor))
	return Spectrum(spectrum.energy_bin_edges[0::factor],
	                spectrum.values[indices].sum(axis=1),
	                (spectrum.errors[indices]**2).sum(axis=1)**(1/2))


def choose_limits(spectrum: "Spectrum", x_label: str, y_label: str) -> tuple[float, float]:
	fig = plt.figure("selection", figsize=FIGURE_SIZE)
	plot_bars(spectrum, x_label, y_label)
	plt.title("click to select the lower and upper bounds of the peak, then close this plot")
	plt.tight_layout()
	lines = [plt.plot([], [], "k--")[0], plt.plot([], [], "k--")[0]]
	curve, = plt.plot(spectrum.energy_bin_edges, np.zeros_like(spectrum.energy_bin_edges),
	                  HIGHLIGHT_COLOR, zorder=1)
	curve.set_visible(False)

	limits: list[float] = []

	def on_click(event: MouseEvent):
		nonlocal limits
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
			if len(limits) == 2:
				gaussian = fit_gaussian(spectrum, min(limits), max(limits))
				curve.set_ydata(gaussian_function(spectrum.energy_bin_edges, gaussian))
				curve.set_visible(True)
			else:
				curve.set_visible(False)
	fig.canvas.mpl_connect('button_press_event', on_click)

	while plt.fignum_exists("selection"):
		plt.pause(.1)
	if len(limits) != 2:
		raise ValueError("you didn't specify both limits.")

	# once the user is done, arrange the results
	return min(limits), max(limits)


def fit_gaussian(spectrum: "Spectrum", left: float, right: float) -> "Gaussian":
	energy_bin_centers = (spectrum.energy_bin_edges[0:-1] + spectrum.energy_bin_edges[1:])/2
	in_limits = (energy_bin_centers >= left) & (energy_bin_centers <= right)
	raw_total = np.sum(spectrum.values*spectrum.energy_bin_widths, where=in_limits)
	popt, pcov = optimize.curve_fit(f=gaussian_function,
	                                xdata=energy_bin_centers[in_limits],
	                                ydata=spectrum.values[in_limits],
	                                sigma=spectrum.errors[in_limits],
	                                p0=[raw_total, (left + right)/2, (right - left)/2])
	values = []
	for i in range(len(popt)):
		values.append(Quantity(popt[i], sqrt(pcov[i, i])))
	return Gaussian(values[0], values[1], values[2])


def gaussian_function(x: NDArray[float],
                      *args: Union["Gaussian", float]) -> NDArray[float]:
	if len(args) == 1:
		N, μ, σ = args[0].total.value, args[0].mean.value, args[0].sigma.value
	elif len(args) == 3:
		N, μ, σ = args
	else:
		raise TypeError("gaussian_function() must be passed a Gaussian object or three floats.")
	return N/sqrt(2*pi*σ**2)*np.exp(-(x - μ)**2/(2*σ**2))


def plot_bars(spectrum: "Spectrum", x_label: str, y_label: str,
              color: str = "k", label: Optional[str] = None) -> None:
	x = np.repeat(spectrum.energy_bin_edges, 2)[1:-1]
	y = np.repeat(spectrum.values, 2)
	plt.plot(x, y, color, linewidth=0.7, zorder=2, label=label)
	energy_bin_centers = (spectrum.energy_bin_edges[:-1] + spectrum.energy_bin_edges[1:])/2
	plt.errorbar(x=energy_bin_centers, y=spectrum.values, yerr=spectrum.errors,
	             ecolor="k", elinewidth=0.7, fmt="none")
	plt.xlim(spectrum.energy_bin_edges[0], spectrum.energy_bin_edges[-1])
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid("on")


def annotate_plot(text: str) -> None:
	text = plt.text(.99, .98, text, zorder=3,
	                ha='right', va='top', transform=plt.gca().transAxes)
	text.set_bbox(dict(facecolor='w', alpha=0.5, edgecolor="none"))


class InvalidFileError(Exception):
	pass


class Spectrum:
	def __init__(self, energy_bin_edges: NDArray[float], values: NDArray[float], errors: NDArray[float]):
		assert energy_bin_edges.ndim == 1 and values.ndim == 1 and errors.ndim == 1
		assert energy_bin_edges.size - 1 == values.size and values.size == errors.size
		self.energy_bin_edges = energy_bin_edges
		self.values = values
		self.errors = errors
		self.energy_bin_widths = energy_bin_edges[1:] - energy_bin_edges[0:-1]


class Quantity:
	def __init__(self, value: float, error: float):
		self.value = value
		self.error = error

	def __mul__(self, other: float) -> "Quantity":
		return Quantity(self.value*other, self.error*other)

	def __format__(self, format_spec: str) -> str:
		# for "e", manage the exponent manually to make the value and error match
		if "e" in format_spec:
			exponent = floor(log10(self.value))
			new_format_spec = format_spec.replace("e", "f")
			return f"{format(self.value/10**exponent, new_format_spec)}e{exponent:+03d} ± " \
			       f"{format(self.error/10**exponent, new_format_spec)}e{exponent:+03d}"
		# otherwise just delegate tothe built-in format for each of the value and error
		else:
			return f"{format(self.value, format_spec)} ± " \
			       f"{format(self.error, format_spec)}"


class Gaussian:
	def __init__(self, total: Quantity, mean: Quantity, sigma: Quantity):
		self.total = total
		self.mean = mean
		self.sigma = sigma


if __name__ == "__main__":
	if len(sys.argv) == 2:
		directory = sys.argv[1:]
		if not os.path.isfile("arguments.txt"):
			with open("arguments.txt", "w") as file:
				file.write(f"directory={directory}\n")
	elif len(sys.argv) == 1:
		directory = None
		try:
			with open("arguments.txt", "r") as file:
				for line in file:
					if "=" in line:
						key = line[:line.index("=")].lower()
						value = line[line.index("=") + 1:].lower().strip()
						if "directory" in key or "path" in key:
							directory = value
		except IOError:
			raise ValueError("You must pass the directory to analyze as an argument.")
		if directory is None:
			raise ValueError("The `arguments.txt` file existed but didn't contain a `directory` argument.")
	else:
		raise ValueError("You must run this script with exactly one command line arguments – "
		                 "the directory to analyze")

	main(directory)
