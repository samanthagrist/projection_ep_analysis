**These scripts are designed to analyze .tif stacks of fluorescence microscopy monitoring images during cell lysis.  Contains:**

1. LysisAnalysisWrapper.m: wrapper script to analyze a set of lysis monitoring stacks
2. Initialization.xlsx: initialization file containing information about each lysis monitoring image set to be analyzed (read in by the wrapper script)
3. A set of various other functions called by the wrapper script during analysis

**Example data that were analyzed by these scripts are available: https://doi.org/10.6078/D1B13V**


**This collection of analysis scripts requires externally-developed code:**

1. The particle-tracking library made available by Blair and Dufresne (based upon the IDL particle tracking software by Grier, Crocker, and Weeks).  Blair and Dufresne's MATLAB implementation available at: http://site.physics.georgetown.edu/matlab/ -- this repository contains trackModified.m (called by the analysis scripts), which is a slightly modified version of Blair and Dufresne's track.m

2. distinguishable_colors.m by Tim Holy, shared on the MATLAB file exchange: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors


