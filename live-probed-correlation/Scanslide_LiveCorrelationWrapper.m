%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%**                 LIVE CELL CORRELATION ANALYSIS                      **%
%**                        October 9th, 2019                            **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

% This wrapper script does image processing analysis on (1) a stitched
% image of live-cell imaging (in-well) and (2) a stitched image of probed
% protein bands (top view) to measure the correlation between live cell
% positions and probed band positions in projection EP.  Both images are
% specified in the configuration .csv file.  

%*************************************************************************%
%**                        SET UP THE PATHS                             **%
%*************************************************************************%
close all
clear
clc

% set the root directory where the confocal image analysis scripts and 
% functions live, set up so that it runs on a Windows and MacOS system 
% with different paths
if ispc % Windows path
    dataAnalysisRoot = 'path1';
    plotSpreadRoot = 'path2';
else % MacOS path
    dataAnalysisRoot = 'path3';
    plotSpreadRoot = 'path4';
end
addpath(dataAnalysisRoot, plotSpreadRoot);

% This script requires the Open Microscopy Bio-Formats MATLAB toolbox,
% available: http://www.openmicroscopy.org/bio-formats/downloads/

% add the bioformats image import code to the path
if ispc
    bioformatsPath = 'path5';
else
    bioformatsPath = 'path6';
end
addpath(bioformatsPath);


%*************************************************************************%
%**                     SET UP THE CONFIGURATIONS                       **%
%*************************************************************************%

DELETEFLAG = intmin(); % use this integer value to flag fields for deletion
% Read in channel configuration information
% load in the configuration .csv file
configFilename = 'Initialization.xlsx';
dataOutputFilename = '2019-10-09_CorrelationAnalysis.mat';
configData = readtable(configFilename);

%numSegments = 4; % the number by which to divide the image in x and y to create ROIs (numSegments^2 ROIs)
roiSz = 175; % size in microns of the length and width of the rectangular ROIs
wBkg = 10; % width of the background 'gutter' region surrounding each ROI
wNoiseBkg = 5; % width of the extra background 'gutter' region surrounding each ROI used for calculation of noise in background-subtracted intensity profiles
wBkgXY = 5; % portion of the X-Y ROI to use as background when measuring x and y intensity profiles for XY resolution calcs
maxNumROIs = 100; % maximum number of ROIs to analyze in each image
fudge = 0.54; % sensitivity for adaptive thresholding to find ROIs
minSpotSize = 30; % minimum size to say we have detected a protein spot
minCellSize = 10; % minimum size to say we have detected a cell
borderSz = 200; % don't segment anything this far from the border
rsqThresh = 0.7; % minimum r squared to accept fitted peak
restartAnalysis = 0;

% if we have already saved the configuration file, load it here.
if exist([dataOutputFilename], 'file') && ~restartAnalysis
    load([dataOutputFilename]);
    configData = readtable(configFilename);
    continuePreviousAnalysis = 1;
end

% Read in list of input and output filenames, and titles
channelConfigFilename = 'ChannelInitialization.xlsx';
channelData = readtable(channelConfigFilename);
numChannels = size(channelData, 1);
proteins = channelData{:, 1};
colourR = channelData{:, 2};
colourG = channelData{:,3};
colourB = channelData{:,4};
revColourR = 1-colourR;
revColourG = 1-colourG;
revColourB = 1-colourB;
defaultColours = [colourR,colourG,colourB];
channelScaling = channelData{:, 5};

numStacks = size(configData, 1);

%% 
%*************************************************************************%
%**                     SET UP THE CONFIGURATIONS                       **%
%*************************************************************************%
% Loop through the list of stacks
for stackIdx = 1:numStacks
%for stackIdx = 11:13
    disp(['Image ', num2str(stackIdx), ' setup'])
    % Set up the filenames
    StackResults(stackIdx).liveFilename = cell2mat(configData{stackIdx, 1});
    StackResults(stackIdx).probedFilename = cell2mat(configData{stackIdx, 2});
    StackResults(stackIdx).outputFilenameBase = cell2mat(configData{stackIdx, 3}); 
    StackResults(stackIdx).liveScale = (configData{stackIdx, 4});
    StackResults(stackIdx).probedScale = (configData{stackIdx, 5}); 
    StackResults(stackIdx).flippedGel = (configData{stackIdx, 6});
    StackResults(stackIdx).roiSz = roiSz;
    StackResults(stackIdx).wBkgXY = wBkgXY;
    StackResults(stackIdx).thresholdFudge = fudge;
    StackResults(stackIdx).rsqThresh = rsqThresh;
    
    save(dataOutputFilename, 'StackResults');
end

%% 
%*************************************************************************%
%**          READ IN AND FIND EXTENTS OF THE STITCHED IMAGES            **%
%*************************************************************************%
% Loop through the list of stacks
for stackIdx = 1:numStacks
%for stackIdx = 11:13
    disp(['Image ', num2str(stackIdx), ' read in and analysis'])
    % Read in the pair of images
    StackResults(stackIdx).liveImage = imread(StackResults(stackIdx).liveFilename);
    StackResults(stackIdx).probedImage = imread(StackResults(stackIdx).probedFilename);
    % Plot the pair of images
    live = sgautoscale(StackResults(stackIdx).liveImage).*channelScaling(1);
    probed = sgautoscale(StackResults(stackIdx).probedImage).*channelScaling(2);
    figure(1); imshow(live);
    figure(2); imshow(probed);
    
    % Allow the user to select 3 corners of the gel on each image.
    % Get the well localization information: (1) array top left, (2) array
    % top right, (3) array bottom left
    % do we need to find the extents?
    if isfield(StackResults(stackIdx), 'extentsFound')
        if ~isempty(StackResults(stackIdx).extentsFound)
            if all(StackResults(stackIdx).extentsFound) > 0
                findExtents = 0;
            else
                findExtents = 1;
            end
        else
            findExtents = 1;
        end
    else
        findExtents = 1;
    end       
    
    if findExtents
        button = 'n';
        % get the gel bounds for the live cell image and probed gel image
        while button == 'n'
            % live cell image
            StackResults = getArrayExtents(StackResults, stackIdx, live, 1);
            % probed gel image
            StackResults = getArrayExtents(StackResults, stackIdx, probed, 2);
            figure(1); title('Press n to redo the array bounds.  Press any other key to confirm.')
            [x,y,button] = ginput(1);
            button
        end
        save(dataOutputFilename, 'StackResults')
    end
end
%% 
%*************************************************************************%
%**              READ IN AND ANALYSE THE STITCHED IMAGES                **%
%*************************************************************************%
% Loop through the list of stacks
for stackIdx = 1:numStacks
    % Rotate and crop the live and probed images according to the bounds we
    % have defined.
    % LIVE
    temp = StackResults(stackIdx).liveImage;
    bottomRight = StackResults(stackIdx).arrayPos{1,1}; % selected as bottom right
    bottomLeft = StackResults(stackIdx).arrayPos{1,2}; % selected as bottom left
    topRight = StackResults(stackIdx).arrayPos{1,3}; % selected as top right
    rotAngle = atand((bottomRight(2)-bottomLeft(2))/(bottomRight(1)-bottomLeft(1)));
    centre = [(size(temp, 1)/2), (size(temp, 2)/2)];
    temp = imrotate(temp, rotAngle);
    R=[cosd(rotAngle) -sind(rotAngle); sind(rotAngle) cosd(rotAngle)];
    coords = [bottomRight; bottomLeft; topRight];
    newCentre = [(size(temp, 1)/2), (size(temp, 2)/2)];
    coordswrtCentre = coords-centre;
    rotcoords = coordswrtCentre*R + newCentre;
    
    bottomRight = rotcoords(1, :);
    bottomLeft = rotcoords(2, :);
    topRight = rotcoords(3, :);
%     figure; imshow(sgautoscale(temp).*channelScaling(1));
%     hold on; plot(rotcoords(:,1), rotcoords(:,2), '.r');
    temp2 = temp(topRight(2):bottomRight(2), bottomLeft(1):bottomRight(1));
%     figure; imshow(sgautoscale(temp2).*channelScaling(1));
    StackResults(stackIdx).cropRotateLiveImage = temp2;
    
    % PROBED
    temp = StackResults(stackIdx).probedImage;
    bottomRight = StackResults(stackIdx).arrayPos{2,1}; % selected as bottom right
    bottomLeft = StackResults(stackIdx).arrayPos{2,2}; % selected as bottom left
    topRight = StackResults(stackIdx).arrayPos{2,3}; % selected as top right
    rotAngle = atand((bottomRight(2)-bottomLeft(2))/(bottomRight(1)-bottomLeft(1)));
    centre = [(size(temp, 1)/2), (size(temp, 2)/2)];
    temp = imrotate(temp, rotAngle);
    R=[cosd(rotAngle) -sind(rotAngle); sind(rotAngle) cosd(rotAngle)];
    coords = [bottomRight; bottomLeft; topRight];
    newCentre = [(size(temp, 1)/2), (size(temp, 2)/2)];
    coordswrtCentre = coords-centre;
    rotcoords = coordswrtCentre*R + newCentre;
    
    bottomRight = rotcoords(1, :);
    bottomLeft = rotcoords(2, :);
    topRight = rotcoords(3, :);
%     figure; imshow(sgautoscale(temp).*channelScaling(2));
%     hold on; plot(rotcoords(:,1), rotcoords(:,2), '.r');
    temp2 = temp(topRight(2):bottomRight(2), bottomLeft(1):bottomRight(1));
    if StackResults(stackIdx).flippedGel
        temp2 = temp2';
    end
%     figure; imshow(sgautoscale(temp2).*channelScaling(2));
    StackResults(stackIdx).cropRotateProbedImage = temp2;
    
    % resample the probed image to align with the live image
    liveSz = size(StackResults(stackIdx).cropRotateLiveImage);
    probedSz = size(StackResults(stackIdx).cropRotateProbedImage);   
    F = griddedInterpolant(double(temp2));
    [sx,sy,sz] = size(temp2);
    xq = (0:probedSz(1)/liveSz(1):sx)';
    yq = (0:probedSz(2)/liveSz(2):sy)';
    %zq = (1:sz)';
    vq = uint16(F({xq,yq}));
    StackResults(stackIdx).scaledCropRotateProbedImage = vq;
    
    % make overlay image
    live = double(StackResults(stackIdx).cropRotateLiveImage);
    probed = double(vq);
    clear('vq');
    probedSz = size(probed);
    overlaySz = [min([liveSz(1), probedSz(1)]), min([liveSz(2), probedSz(2)])];
    overlayFig = uint8(zeros(overlaySz(1), overlaySz(2), 3));
    overlayFig(:,:,1) = overlayFig(:,:,1) + revColourR(1).*sgautoscale(live(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(1) + revColourR(2).*sgautoscale(probed(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(2);
    overlayFig(:,:,2) = overlayFig(:,:,2) + revColourG(1).*sgautoscale(live(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(1) + revColourG(2).*sgautoscale(probed(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(2);
    overlayFig(:,:,3) = overlayFig(:,:,3) + revColourB(1).*sgautoscale(live(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(1) + revColourB(2).*sgautoscale(probed(1:overlaySz(1), 1:overlaySz(2))).*channelScaling(2);
    figure; imshow(overlayFig);
    imwrite(overlayFig, [StackResults(stackIdx).outputFilenameBase, 'OverlayImage'],'tif');
    imwrite(overlayFig, [StackResults(stackIdx).outputFilenameBase, 'OverlayImage'],'png');
    overlayFigRev = imcomplement(overlayFig);
    hOverlay = figure; imshow(overlayFigRev);
    imwrite(overlayFigRev, [StackResults(stackIdx).outputFilenameBase, 'OverlayImageInvert'],'tif');
    imwrite(overlayFigRev, [StackResults(stackIdx).outputFilenameBase, 'OverlayImageInvert'],'png');
    
    save(dataOutputFilename, 'StackResults');
       
    % Segment single cells from the live cell image 
    tmp = segmentImage(StackResults(stackIdx), live, 0.35, minCellSize, borderSz);
    StackResults(stackIdx).liveCentroids = tmp.centroids;
    StackResults(stackIdx).liveROILocX = tmp.ROILocX;
    StackResults(stackIdx).liveROILocY = tmp.ROILocY;
    StackResults(stackIdx).numLiveROIs = tmp.numROIs;
    
    % Segment protein spots from the probed gel image 
    tmp = segmentImage(StackResults(stackIdx), probed, 0.55, minSpotSize, borderSz);
    StackResults(stackIdx).probedCentroids = tmp.centroids;
    StackResults(stackIdx).probedROILocX = tmp.ROILocX;
    StackResults(stackIdx).probedROILocY = tmp.ROILocY;
    StackResults(stackIdx).numProbedROIs = tmp.numROIs;
 
    
    % Plot spatial map of live cell and probed band positions
    h = figure; plot(StackResults(stackIdx).liveROILocX, overlaySz(2)-StackResults(stackIdx).liveROILocY, '.', 'Color', [colourR(1), colourG(1), colourB(1)]);
    hold on; plot(StackResults(stackIdx).probedROILocX, overlaySz(2)-StackResults(stackIdx).probedROILocY, '.', 'Color', [colourR(2), colourG(2), colourB(2)]);
    saveas(h, [StackResults(stackIdx).outputFilenameBase, 'segCoords'], 'fig');
    saveas(h, [StackResults(stackIdx).outputFilenameBase, 'segCoords'], 'eps');
    
    % Identify probed bands that are within some fudge region of live cells. 
    StackResults(stackIdx).matchedLiveSegs = ismembertol([StackResults(stackIdx).liveROILocX, StackResults(stackIdx).liveROILocY], ...
        [StackResults(stackIdx).probedROILocX, StackResults(stackIdx).probedROILocY], StackResults(stackIdx).roiSz*StackResults(stackIdx).liveScale, ...
        'ByRows', true, 'DataScale', 1);
    StackResults(stackIdx).matchedProbedSegs = ismembertol([StackResults(stackIdx).probedROILocX, StackResults(stackIdx).probedROILocY], ...
        [StackResults(stackIdx).liveROILocX, StackResults(stackIdx).liveROILocY], StackResults(stackIdx).roiSz*StackResults(stackIdx).liveScale, ...
        'ByRows', true, 'DataScale', 1);
    StackResults(stackIdx).numMatchedLiveSegs = nnz(StackResults(stackIdx).matchedLiveSegs);
    StackResults(stackIdx).numMatchedProbeSegs = nnz(StackResults(stackIdx).matchedProbedSegs);
    
    
    % What proportion of bands correlate with live cells?
    StackResults(stackIdx).fracMatchedProbeSegs = StackResults(stackIdx).numMatchedProbeSegs/length(StackResults(stackIdx).probedROILocX);
    
    % What proportion of cells correlate with detectable bands?
    StackResults(stackIdx).fracMatchedLiveSegs = StackResults(stackIdx).numMatchedLiveSegs/length(StackResults(stackIdx).liveROILocX);
    
    % Mark coincident spots on overlay image
    figure(hOverlay);
    hold on; plot(StackResults(stackIdx).liveROILocX(StackResults(stackIdx).matchedLiveSegs~=0), ...
        StackResults(stackIdx).liveROILocY(StackResults(stackIdx).matchedLiveSegs~=0), 'ko', 'MarkerSize', 10);
    saveas(hOverlay, [StackResults(stackIdx).outputFilenameBase, 'matchCoordsOverlay'], 'fig');
    saveas(hOverlay, [StackResults(stackIdx).outputFilenameBase, 'matchCoordsOverlay'], 'png');
    
    save(dataOutputFilename, 'StackResults');
    
end










