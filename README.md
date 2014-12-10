stellarstructure
================

Stellar Structure Model for UCSC Astro 220A
-------

Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   

Copyright (c) 2014
Licensed under the MIT License

-------

##Main Objectives
  - Given boundary conditions, integrate from the surface and center of a star to an arbitrary fitting point
    - see: integrate.py
  - Given initial boundary conditions, iteratively fit new boundary conditions that minimize the error at the fitting point between the inward and outward integrations
  - Plot results
    - see: plotmodel.py
  - Print and plot code diagnostics
    - functions incorporated into each module
  - Summarize the work in a written report
    - see: report/Gentry_StellarModel_Report
 

##Requires
  - Python2
    - Numpy
    - Scipy
    - matplotlib
  - At least 2 available CPUs
    - paralellized code within shootf can easily be serialized if necessary
  
  
##Optional Packages
  - LaTeX
