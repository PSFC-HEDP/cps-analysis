# cps-analysis

 Simple command-line tool to extract mean energies and yields from CPS scans.
 
## How to install

 It's just a couple of python scripts and some CSVs; all you need to do is download the repo.

 There are some dependencies.
 Most of them are on PyPI and can be installed with one command if you don't already have them:
 ~~~bash
 pip install matplotlib numpy pandas scipy xarray
 ~~~

 This program also requires Peter Heuer’s CR39py libary.
 It’s not on PyPI last I checked, but you can pip-install it using
 ~~~bash
 pip install git+https://github.com/pheuer/CR39py.git
 ~~~

## How to use

### Infer spectra

 To run the basic spectrum inference, simply navigate to the root directory in a command line and call
 ~~~bash
 python analyze_cps.py {cps1-finger} {cps2-finger} {particle} {path/to/files/}
 ~~~
 where `{cps1-finger}` stands for the CPS1 configuration (for example “b5”),
 `{cps2-finger}` stands for the CPS2 configuration (for example “d9w”),
 `{particle}` stands for the particle being measured or the A/Z^2 of the particle being measured (for example “d” or “2.0”)
 and `{path/to/files/}` stands for the directory, absolute or relative, that contains the `.cpsa` files.
 If you only have files for one CPS, you may set the finger of the other one to “none”.

 You will first be prompted to click on the plot to highlight the background region.
 This will generally be somewhere in the empty space above the indicated data region.
 You may specify any two diagonal corners of the rectangle in either order.
 At any time, you may right-click to delete the most recently placed corner.
 When both points are set, the rectangle will be displayed for confirmation.
 Close the plot at this point to continue.

 You will then be prompted to click on the plot to select the diameter limits.
 Each time you click, it will place a new point on the plot.
 To use diagonal diameter cuts, you’ll want to click multiple times from left to right,
 drawing a line either above or below the signal region.
 Once you get all the way to the right, start clicking from the left again to place the limit on the other side.
 Pay attention to how the lines appear to make sure it’s understanding you correctly.
 As before, you can right-click at any point to undo an action.
 Close the plot when the diameter cuts are satisfactory.
 
 It will repeat this for every scan file in the given directory,
 outputing a CSV containing the spectrum for each one.

### Infer yield and mean energy

 To extract yields, mean energies, and plots from the analyzed spectra, call
 ~~~bash
 python analyze_spectra.py {path/to/files}
 ~~~
 where `{path/to/files/}` stands for the directory, absolute or relative, that contains the `.cpsa` files.
 If the directory is ommitted, it will reuse the same directory you specified for the previous step.

 You will be prompted to click on the plot to select the upper and lower limits for the peak fit.
 As before, you can right-click at any time to undo an action.
 Close the plot when the peak limits are satisfactory.
 
 It will repeat this for every spectrum file in the given directory,
 outputing a plot of each spectrum as a PNG,
 with the yield, mean energy, and width written on the side.
 Lmk if you want the output as EPS or something for whatever reason;
 it's not that difficult a thing to add.

 This program isn’t as powerful as Fredrick’s AnalyzeCR39,
 so if you have a use case that doesn’t seem to be covered by the choices you’re given,
 you may need to go use that.
 But I hope that this will cover most use cases in a more convenient format.

## Limitations

 The main limitation is that this script doesn’t make it very easy for the user to go back and forth on their limits.
 Once you set a cut or highlight a region, it will automatically move onto the next question/piece.
 I hope that the fact that this program makes the data much easier to visualize than AnalyzeCR39
 mitigates this shortcoming. 

 This analysis also currently assumes a slit width of 2mm.
 You will need to modify the script or ask me to add an option if you want to use it for 1mm.
 
 The broadening due to the finite slit width is not accounted for when the peak width is calculated,
 so be wary of that if you care about the peak width.

 I don't currently allow the user to set the data region,
 so if part of the scan is bad there's no way to remove it.
 I probably could make it so clicking after you set the background would do it...
 but I haven't.

 Lastly, there is no provision for ramped background subtraction,
 so users will have to make do with flat background for now.

## To do

 In the future, I would really like to have this also automatically compute ρR with error bars
 using the fit mean energy and particle species.
 But I haven't done that yet.
