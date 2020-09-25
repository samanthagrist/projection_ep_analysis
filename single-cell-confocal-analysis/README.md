****These scripts are designed to analyze a set of confocal microscopy images of single-cell projection electrophoresis separations.  Here, the user manually selects the ROIs corresponding to separation lanes in each image.  Contains:**

1. ManualROI_StackAnalysis_autoBkg_MigrationAnalysis_makeXZPlots.m: wrapper script to analyze a set of confocal image stacks
2. Initialization.xlsx: initialization file containing information about the confocal images to be analyzed (read in by the wrapper script)
3. ChannelInitialization.xlsx: initialization file containing information about the different colour channels in each image (read in by the wrapper script)
3. A set of various other functions called by the wrapper script during analysis

**Example data that were analyzed by these scripts are available: https://doi.org/10.6078/D1B13V**


**This collection of analysis scripts requires externally-developed code:**

1. Functions from Bio-Formats (The Open Microscopy Environment bfmatlab): https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html)

2. PlotSpread.m by Jonas, shared on the MATLAB file exchange.  Available: https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot

3. distinguishable_colors.m by Tim Holy, shared on the MATLAB file exchange: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
