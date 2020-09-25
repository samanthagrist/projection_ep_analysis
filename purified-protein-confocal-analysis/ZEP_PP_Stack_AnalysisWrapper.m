%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%**                ZEP STACK ANALYSIS WRAPPER SCRIPT                    **%
%**                         MARCH 1ST, 2018                             **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

% This wrapper script does image processing analysis on a series of Z stacks 
% specified in the configuration .csv file.  Functionality to track
% purified protein ROIs through Z, extract the migration distance and
% peak widths in x, y, and z, and plot.


close all;
clear;

%*************************************************************************%
%**                        SET UP THE PATHS                             **%
%*************************************************************************%

% set the root directory where the confocal image analysis scripts and 
% functions live, set up so that it runs on a Windows and MacOS system 
% with different paths
if ispc % Windows path
    dataAnalysisRoot = 'path1';
else % MacOS path
    dataAnalysisRoot = 'path2';
end
addpath(dataAnalysisRoot);

% This script requires the tracking library written by Daniel Blair and
% Eric Dufresne, available: http://site.physics.georgetown.edu/matlab/index.html

% add the Blair/Dufresne tracking code to the path
addpath([dataAnalysisRoot,'/tracking/'])

% This script also requires the Open Microscopy Bio-Formats MATLAB toolbox,
% available: http://www.openmicroscopy.org/bio-formats/downloads/

% add the bioformats image import code to the path
addpath([dataAnalysisRoot,'/bfmatlab/'])

%*************************************************************************%
%**         INITIALIZATION - SET PARAMETERS FOR THESE DATASETS          **%
%*************************************************************************%

% load in the configuration .csv file
configFilename = '2018-07-23_Initialization.xlsx';
dataOutputFilename = '2018-07-23_Analysis2_';
configData = readtable(configFilename);
numStacks = size(configData, 1);

% set configuration parameters that hold for all image stacks in this set
ESTSZ = 1; % estimate the area of the protein spot?
IThresh = 1e3; % full-image intensity threshold to decide whether to calculate the segmented size in an image
proteins = [{'BSA'}, {'OVA'}, {'IgG'}]; % proteins to analyze
peaknames = [{'BSA'}, {'BSA`'}; ... % all peaknames to analyze in first colour channel
             {'OVA'}, {'OVA`'}; ... % all peaknames to analyze in second colour channel
             {'IgG'}, {''}];          % all peaknames to analyze in third colour channel
colIdx = 'rgk'; % colour order for plotting
halfROISz = 150; % size in pixels to offset from centroid to determine ROI
bkgDilation = 100; % size in pixels of the distance from the thresholded region to use as background (how far away from segmented features)
restartAnalysis = 0; % set to 1 to restart (not load previous data from the output file and OVERWRITE)
zSpacing = 5; % spacing of the z slices (in um)
xyScale = 0.34; % x-y image scale (in um/pixel)

%*************************************************************************%
%**                     INITIALIZATION - AUTOMATIC                      **%
%*************************************************************************%
% initialize the number of found intensity peaks to zero and define the
% labels
CombinedResults.IntensityPeakEntries = 0;
CombinedResults.IntensityDataLabels = [{'Gel number'}, {'Gel %T'}, {'Gel %R'}, {'EP time'}, ...
    {'EP field'}, {'Stack index'}, {'ROI index'}, {'Protein channel'}, {'Peak name'}, ...
    {'Peak A'}, {'Peak Mu'}, {'Peak Sigma'}, {'Peak Rsq'}, {'x res at peak'}, {'y res at peak'}];

% set up the configuration parameters for each image stack in this set
for stackIdx = 1:numStacks
    % set parameters from the csv file data
    StackResults(stackIdx).inputFilename = cell2mat(configData{stackIdx, 1});
    StackResults(stackIdx).outputFilenameBase = cell2mat(configData{stackIdx, 2}); 
    
    StackResults(stackIdx).halfROISz = halfROISz; % size in pixels to offset from centroid to determine ROI
    
    % Initialize the elements we add to zero
    StackResults(stackIdx).bkgMean = 0;
    StackResults(stackIdx).bkgStd = 0;
    StackResults(stackIdx).currProtein = 0;
    StackResults(stackIdx).A = 0;
    StackResults(stackIdx).AScaled = 0;
    StackResults(stackIdx).centroids = {0};
    StackResults(stackIdx).imgSz = 0;
    StackResults(stackIdx).trackSetArray = {0};
    StackResults(stackIdx).trackCoordsArray = {[]};
    StackResults(stackIdx).numGoodParticles = 0;
    StackResults(stackIdx).numZ = 0;
    StackResults(stackIdx).intensityResults = 0;
    StackResults(stackIdx).summedIntensity = 0;
    StackResults(stackIdx).zLoc = 0;
    StackResults(stackIdx).xyFitA = 0;
    StackResults(stackIdx).xyFitMu = 0;
    StackResults(stackIdx).xyFitSigma = 0;
    StackResults(stackIdx).xyFitRSq = 0;
    StackResults(stackIdx).avgIntensity = 0;
    StackResults(stackIdx).stdIntensity = 0;
    StackResults(stackIdx).avgWidth = 0;
    StackResults(stackIdx).stdWidth = 0;
    StackResults(stackIdx).intensityFitBounds = 0;
    StackResults(stackIdx).intPkFitA = 0;
    StackResults(stackIdx).intPkFitMu = 0;
    StackResults(stackIdx).intPkFitSigma = 0;
    StackResults(stackIdx).intPkFitRSq = 0;
    StackResults(stackIdx).intPkFitXYRes = 0;
end
 
    % if we have already saved the configuration file, load it here.
    if exist([dataOutputFilename, 'DATA.mat'], 'file') && ~restartAnalysis
        load([dataOutputFilename, 'DATA.mat']);
        configData = readtable(configFilename);
        continuePreviousAnalysis = 1;
    end
    
% set up the configuration parameters for each image stack in this set
for stackIdx = 1:numStacks   
    % set parameters from the csv file data  
    StackResults(stackIdx).inputFilename = cell2mat(configData{stackIdx, 1});
    StackResults(stackIdx).outputFilenameBase = cell2mat(configData{stackIdx, 2}); 
    StackResults(stackIdx).exampleProtein = configData{stackIdx, 3}; % protein channel to use for ROI tracking segmentation
    StackResults(stackIdx).topZ = configData{stackIdx, 4}; % Z position at the top of the gel/bottom of the well
    StackResults(stackIdx).endZ = configData{stackIdx, 12}; % Last Z position to analyze (set to zero to go through the whole stack)
    StackResults(stackIdx).flipStack = configData{stackIdx, 5}; % set to 1 if the wells are on the high end of the image slices (if we started imaging at the non-well side)
    StackResults(stackIdx).epTime = configData{stackIdx, 6}; % electrophoresis time
    StackResults(stackIdx).epField = configData{stackIdx, 7}; % electric field strength (V/cm) for electrophoresis
    StackResults(stackIdx).gelT = configData{stackIdx, 8}; % gel %T
    StackResults(stackIdx).gelR = configData{stackIdx, 9}; % gel %Rhinohide
    StackResults(stackIdx).thresholdFudge = configData{stackIdx, 10}; % multiplicative factor for Otsu's grey threshold for segmentation
    StackResults(stackIdx).gelNumber = configData{stackIdx, 11}; % the number of the gel used for this condition
    StackResults(stackIdx).tDiff = configData{stackIdx, 13}; % the diffusion time for this experiment
    StackResults(stackIdx).currZ = 0;
  
    StackResults(stackIdx).imgReader = bfGetReader(StackResults(stackIdx).inputFilename); % use the bioformats toolbox to create the image reader for this stack
    reader = StackResults(stackIdx).imgReader;
    reader.setSeries(0);
    planes = reader.getImageCount;
    StackResults(stackIdx).numZ = planes/length(proteins);
    if isfield(StackResults(stackIdx), 'endZ')
        if (StackResults(stackIdx).endZ ~= 0)
            StackResults(stackIdx).numZ = StackResults(stackIdx).endZ;
        end
    end
    reader.close();

    % Set common parameters defined above
    StackResults(stackIdx).ESTSZ = ESTSZ; % estimate the area of the protein spot?
    StackResults(stackIdx).IThresh = IThresh; % full-image intensity threshold to decide whether to calculate the segmented size in an image
    StackResults(stackIdx).proteins = proteins; % proteins to analyze    
    StackResults(stackIdx).peaknames = peaknames; % proteins to analyze
    StackResults(stackIdx).colIdx = colIdx; % colour order for plotting
    StackResults(stackIdx).halfROISz = halfROISz; % size in pixels to offset from centroid to determine ROI
    StackResults(stackIdx).bkgDilation = bkgDilation; % size in pixels to offset from centroid to determine ROI
    StackResults(stackIdx).zSpacing = zSpacing; % spacing of the z slices (in um)
    StackResults(stackIdx).xyScale = xyScale; % x-y image scale (in um/pixel)
    StackResults(stackIdx).stackIdx = stackIdx; % the index of this image

    
end


%*************************************************************************%
%**                  EXTRACT INTENSITY INFORMATION                      **%
%*************************************************************************%
for stackIdx = 1:numStacks
    close all;
    StackResults(stackIdx).imgReader = bfGetReader(StackResults(stackIdx).inputFilename); % use the bioformats toolbox to create the image reader for this stack

    % if the background mean and stdev are not yet defined, call the getBkg
    % function.
    stackIdx
    if (StackResults(stackIdx).bkgMean(1) == 0) || (StackResults(stackIdx).bkgStd(1) == 0)
        StackResults(stackIdx) = getbkg(StackResults(stackIdx));
        save([dataOutputFilename, 'DATA.mat'], 'StackResults', 'CombinedResults');
    end
    
    % if we haven't yet tracked the centroids through Z for this stack,
    % call the trackcentroids function
    if StackResults(stackIdx).numGoodParticles == 0
        StackResults(stackIdx) = trackcentroids(StackResults(stackIdx));
        save([dataOutputFilename, 'DATA.mat'], 'StackResults', 'CombinedResults');
    end
    
    % if we haven't yet calculated the summed ROI intensities through Z,
    % call the calcintensities function
    if StackResults(stackIdx).intensityResults == 0
        StackResults(stackIdx) = calcintensities(StackResults(stackIdx));
        save([dataOutputFilename, 'DATA.mat'], 'StackResults', 'CombinedResults');
    end
    
    % plot the calculated intensity profiles
    if StackResults(stackIdx).avgIntensity == 0
        StackResults(stackIdx) = plotintensities(StackResults(stackIdx));
    end
    StackResults(stackIdx).imgReader.close()
end

%% 


%*************************************************************************%
%**               MEASURE AND PLOT MIGRATION INFORMATION                **%
%*************************************************************************%
for stackIdx = 1:numStacks
    % if we haven't yet performed the migration analysis on this dataset,
    % call the measuremigration function
    if numel(StackResults(stackIdx).intPkFitMu) < 2  
        %CombinedResults.IntensityPeaks = {};
        [StackResults(stackIdx), CombinedResults] = measuremigration(StackResults(stackIdx), CombinedResults); 
        save([dataOutputFilename, 'DATA.mat'], 'StackResults', 'CombinedResults');
    end
    
    
end

% plot the migration information
[StackResults(stackIdx), CombinedResults] = plotmigration(StackResults(stackIdx), CombinedResults);






    
    