protoDUNE Cosmic Ray Muons
==========================

This repo has files for generating cosmic ray muons. The generators simply
sample a parameterized distribution. Events our output in HEPEVT text format.
The first ROOT script version of the software is documented in LBNE doc-db
9647, and is by Clay Barton.

SamplingProgram.C is the main script for the first version. As far as I can
tell, it samples theta from cos(theta)^2 and then samples from the E
distribution.

The second version is by Justin Hugon, and the main script is mycosmics.py.
mycosmics.py samples from the 2D distribution w.r.t. theta and E.
mycosmics_differentialplots.py plots the differential function used in
mycosmics.py and directly integrates it. 

There is also a converter from HEPEVT to ROOT tree format program:
makeRootTree.py.  It requires the rootpy python package.

Other references:
-----------------

Mengyun, Guan, (2011). Muon Simulation at the Daya Bay Site, Lawrence
Berkeley National Laboratory, LBNL Paper LBNL-4262E. Retrieved from:
http://escholarship.org/uc/item/6jm8g76d

Chirkin, Dmitry (2004). Fluxes of Atmospherice Leptons at 600 GeV - 60 TeV,
hep-ph 0407078.

Gaisser, Thomas and Stanev, Todor (2004). Cosmic Rays, Review Of Particle
Physics, Physics Letter B 592.
