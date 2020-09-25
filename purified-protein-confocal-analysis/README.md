**These scripts are designed to analyze confocal images of purified protein projection electrophoresis separations.  Contains:**

1. ZEP_PP_Stack_AnalysisWrapper.m: wrapper script to analyze a set of purified protein confocal images
2. 2018-07-23_Initialization.xlsx: initialization file containing information about each confocal image to be analyzed (read in by the wrapper script)
3. A set of various other functions called by the wrapper script during analysis

**Example data that were analyzed by these scripts are available: https://doi.org/10.6078/D1B13V**


**This collection of analysis scripts requires externally-developed code:**

1. Functions from Bio-Formats (The Open Microscopy Environment bfmatlab): https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html)

2. The particle-tracking library made available by Blair and Dufresne (based upon the IDL particle tracking software by Grier, Crocker, and Weeks).  Blair and Dufresne's MATLAB implementation available at: http://site.physics.georgetown.edu/matlab/ -- this repository contains trackModified.m (called by the analysis scripts), which is a slightly modified version of Blair and Dufresne's track.m
