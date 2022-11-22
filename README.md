# cps-analysis

 Simple compand-line utility to extract charged particle spectra from CPS scans.

 To run the analysis, simply navigate to the root directory in the command line and call
 ~~~bash
    python analyze_cps.py {cps1-finger} {cps2-finger} path/to/files/
 ~~~
 where `{cps1-finger}` stands for the CPS1 configuration (for example “b8w”),
 `{cps2-finger}` stands for the CPS2 configuration (for example “d9”),
 and `path/to/files/` stands for the directory, absolute or relative, that contains the `.cpsa` files.
 If one of the CPS was not run, you may set its configuration to “none”.
 Follow the instructions on screen to choose diameter cuts and such.

 This program isn’t as powerful as Fredrick’s AnalyzeCR39,
 so if you have a use case that doesn’t seem to be covered by the choices you’re given,
 you may need to go use that.
 But I hope that this will cover most use cases in a more convenient format.

## Dependencies

 This program requires Peter Heuer’s CR39py libary.
 I’ll drop in some instructions for how to install that later.
