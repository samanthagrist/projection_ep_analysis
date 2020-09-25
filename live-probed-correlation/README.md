**These scripts are designed to measure the correlation between live cell images and probed wide field fluorescence images of projection EP gels after separation and immunoprobing.  Contains:**

1. Scanslide_LiveCorrelationWrapper.m: wrapper script to analyze a pair of live/probed images
2. Initialization.xlsx: initialization file containing information about the two stitched fluorescence images to be analyzed (read in by the wrapper script)
3. ChannelInitialization.xlsx: initialization file containing false-colour information about the fluorescence contained in each image (read in by the wrapper script)
3. A set of various other functions called by the wrapper script during analysis

**Example data that were analyzed by these scripts are available: https://doi.org/10.6078/D1B13V**


**This collection of analysis scripts requires externally-developed code:**

1. Functions from Bio-Formats (The Open Microscopy Environment bfmatlab): https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html)

2. PlotSpread.m by Jonas, from the MATLAB file exchange.  Available: https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
