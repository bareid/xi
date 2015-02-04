xi
==

This is a work-in-progress module to compute 2-point correlation functions on both survey data and in periodic box simulations.  Input bin files allow a variety of two-point correlation functions to be computed (projected, multipoles, rp-pi diagrams, etc.)

If you want to reproduce the binning scheme (``xiellhat'') used in the Reid et al. 2014 analysis of small-scale CMASS clustering:

1) Run xi with binfile xibinfiles/bin1xiellsmallscale.txt
2) Input the result to xiell.xiellfromDR(fbase=YOURFBASE,binfile=bin1fineMU.txt,rperpcut=0.533655) from my LSSanalysis package.  Add periodicopt=1 if you ran xi with the sims option (periodic boundary conditions).  The input bin1fineMU.txt file is available with the LSSanalysis package.

The rebinning option should eventually be incorporated into xi itself, but for now, this scheme works.

