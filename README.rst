protoDUNE Cosmic Ray Muons
==========================

This repo has files for generating cosmic ray muons. The generator simply
uses accept/reject sampling. Events are output in HEPEVT text format.
The first ROOT script version of the software is documented in LBNE doc-db
9647, and is by Clay Barton. The second version is by Justin Hugon, and the
main script is mycosmics.py.

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
