# cps-analysis

 Simple compand-line utility to extract charged particle spectra from CPS scans.

## How to use

 To run the basic, simply navigate to the root directory in the command line and call
 ~~~bash
 python analyze_cps.py {cps1-finger} {cps2-finger} {particle} {path/to/files/}
 ~~~
 where `{cps1-finger}` stands for the CPS1 configuration (for example “b8w”),
 `{cps2-finger}` stands for the CPS2 configuration (for example “d9”),
 `{particle}` stands for the particle being measured or the A/Z^2 of the particle being measured (for example “d” or “2”)
 and `{path/to/files/}` stands for the directory, absolute or relative, that contains the `.cpsa` files.
 If one of the CPS was not run, you may set its configuration to “none”.

 You will first be prompted to click on the plot to highlight the data region.
 This will generally be somewhere in the space above the indicated data region.
 You may specify any two diagonal corners of the rectangle in either order.
 At any time, you may right-click to delete the most recently placed point.
 When both points are set, the rectangle will be displayed for confirmation.
 Close the plot at this point to continue.

 You will be prompted to click on the plot to select the diameter limits.
 Each time you click, it will place a new point on the plot.
 To use diagonal diameter cuts, you’ll want to click from left to right,
 drawing a line either above or below the signal region.
 Once you get to the right end, start clicking from the left again to place the limit on the other side.
 Pay attention to how the lines appear to make sure it’s understanding you correctly.
 As before, you can right-click at any point to undo an action.
 Close the plot when the diameter cuts are satisfactory.

 To extract yields, means, and plots from analyzed spectra,
 from the root directory, call the follow-up script
 ~~~bash
 python analyze_spectra.py {path/to/files}
 ~~~
 where `{path/to/files/}` stands for the directory, absolute or relative, that contains the `.cpsa` files.
 If the directory is ommitted, the same directory used for the previous script will be assumed.

 This will prompt you to click on the plot to select limits for the peak.
 As before, you can right-click at any time to undo an action.
 Close the plot when the peak limits are satisfactory.

 This program isn’t as powerful as Fredrick’s AnalyzeCR39,
 so if you have a use case that doesn’t seem to be covered by the choices you’re given,
 you may need to go use that.
 But I hope that this will cover most use cases in a more convenient format.

## Limitations

 The main limitation is that this script doesn’t make it very easy for the user to go back and forth on their limits.
 Once you set a cut or highlight a region, it will automatically move onto the next question/piece.
 I hope that the fact that this program makes the data much easier to visualize than AnalyzeCR39
 mitigates this shortcoming somewhat. 

 This analysis also currently assumes a slit width of 2mm.
 You will need to modify the script if you want to use it for 1mm.

 I don't currently allow the user to set the data region,
 so if part of the scan is bad there's no way to remove it.
 I probably could make it so clicking after you set the background would do it...
 but I haven't.

 Lastly, there is no provision for ramped background subtraction,
 so users will have to make do with flat background for now.

## Dependencies

 This program requires Peter Heuer’s CR39py libary.
 It’s not on PyPI last I checked, but you can pip-install it using
 ~~~bash
 pip install git+https://github.com/pheuer/CR39py.git
 ~~~
