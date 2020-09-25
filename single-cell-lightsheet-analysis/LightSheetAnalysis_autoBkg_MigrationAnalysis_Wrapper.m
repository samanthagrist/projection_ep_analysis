%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%**                   LIGHT SHEET ANALYSIS SCRIPT                       **%
%**                      September 12th, 2019                           **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

% This wrapper script does image processing analysis on a tiled Z stack 
% specified in the configuration .csv file.  Functionality to track
% protein ROIs through Z, extract the migration distance and
% peak widths in x, y, and z, and plot.

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
dataOutputFilename = '2019-10-04_Analysis_Gel5_Tiled_CorrectedTopZ.mat';
configData = readtable(configFilename);
% check that we only have a single stack set up to be analyzed, and send a
% warning if more were entered in the config file
numStacks = size(configData, 1);
if numStacks ~= 1
    warning('Light sheet analysis script only analyzes a single .czi file.  All additional files in the configuration list will be ignored.')
end
%numSegments = 4; % the number by which to divide the image in x and y to create ROIs (numSegments^2 ROIs)
roiSz = 180; % size in pixels of the length and width of the rectangular ROIs
wBkg = 10; % width of the background 'gutter' region surrounding each ROI
wNoiseBkg = 5; % width of the extra background 'gutter' region surrounding each ROI used for calculation of noise in background-subtracted intensity profiles
wBkgXY = 5; % portion of the X-Y ROI to use as background when measuring x and y intensity profiles for XY resolution calcs
maxNumROIs = 100; % maximum number of ROIs to analyze in each image
fudge = 0.56; % sensitivity for adaptive thresholding to find ROIs
minSpotSize = 200; % minimum size to say we have detected a protein spot
rsqThresh = 0.7; % minimum r squared to accept fitted peak
restartAnalysis = 1;


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
defaultColours = [colourR,colourG,colourB];
channelScaling = channelData{:, 5};

% Set up the reader and extract the number of tiles in the .czi stack
stackIdx = 1; % only look at the first entry of the config data file for light sheet
inputFilename = cell2mat(configData{stackIdx, 1});
reader = bfGetReader(inputFilename); % use the bioformats toolbox to create the image reader for this stack
data = reader.getMetadataStore(); % get metadata
numTiles = reader.getSeriesCount; % get the number of series (tiles) in this image set

%% 


%*************************************************************************%
%**                     SET UP THE CONFIGURATIONS                       **%
%*************************************************************************%
% Loop through the list of stacks
for tileIdx = 1:numTiles
%for tileIdx = 11:13
    disp(['Tile ', num2str(tileIdx)])
    % Set up the filenames
    StackResults(tileIdx).proteins = proteins;
    StackResults(tileIdx).defaultColours = defaultColours;
    StackResults(tileIdx).inputFilename = cell2mat(configData{stackIdx, 1});
    StackResults(tileIdx).outputFilenameBase = cell2mat(configData{stackIdx, 2}); 
    StackResults(tileIdx).title = cell2mat(configData{stackIdx, 3});
    StackResults(tileIdx).flipStack = configData{stackIdx, 4};
    StackResults(tileIdx).topZ = configData{stackIdx, 5};
    StackResults(tileIdx).segChannel = configData{stackIdx, 6}; % the channel to use for the segmentation.  Set to zero to use sum of all
    StackResults(tileIdx).roiSz = roiSz;
    StackResults(tileIdx).wBkgXY = wBkgXY;
    StackResults(tileIdx).thresholdFudge = fudge;
    StackResults(tileIdx).rsqThresh = rsqThresh;
    
    numPeaks = size(StackResults(tileIdx).proteins, 2);

    % Get the XY and Z resolution (value and unit scale factor)
    StackResults(tileIdx).voxelSizeXYdefaultValue = data.getPixelsPhysicalSizeX(tileIdx - 1).value();
    StackResults(tileIdx).voxelSizeZdefaultValue = data.getPixelsPhysicalSizeZ(tileIdx - 1).value();
    StackResults(tileIdx).voxelSizeXYunitScale = data.getPixelsPhysicalSizeX(tileIdx - 1).unit().getScaleFactor; % get the unit scale factor to get to [m]
    StackResults(tileIdx).voxelSizeZunitScale = data.getPixelsPhysicalSizeZ(tileIdx - 1).unit().getScaleFactor;
    
    % Get the stage positions
    planeIndex = 0; % start with the first image in the stack to get some example data (value and unit scale factor)
    StackResults(tileIdx).stageX = data.getPlanePositionX(tileIdx - 1, planeIndex).value().doubleValue; % get the image plane stage X location
    StackResults(tileIdx).stageY = data.getPlanePositionY(tileIdx - 1, planeIndex).value().doubleValue; % get the image plane stage Y location
    StackResults(tileIdx).stageXUnitScale = data.getPlanePositionX(tileIdx - 1, planeIndex).unit().getScaleFactor; % get the unit scale factor to get to [m]
    StackResults(tileIdx).stageYUnitScale = data.getPlanePositionY(tileIdx - 1, planeIndex).unit().getScaleFactor;
    relativeScaleX = StackResults(tileIdx).voxelSizeXYunitScale.doubleValue/StackResults(tileIdx).stageXUnitScale.doubleValue; % find the relative scale between the pixel resolution and stage units
    relativeScaleY = StackResults(tileIdx).voxelSizeXYunitScale.doubleValue/StackResults(tileIdx).stageYUnitScale.doubleValue; % find the relative scale between the pixel resolution and stage units
    % Transform the stage coordinates to be in the same units as the pixel
    % resolution if they are not already
    StackResults(tileIdx).stageX = StackResults(tileIdx).stageX/relativeScaleX;
    StackResults(tileIdx).stageY = StackResults(tileIdx).stageY/relativeScaleY;

    % Check that we have the right number of colour channels
    if data.getChannelCount(tileIdx-1) ~= numChannels
        warning('The number of defined colour channels is not equal to the image colour channels');
        numChannels = min([data.getChannelCount(tileIdx-1), numChannels]);
    end
    StackResults(tileIdx).numChannels = numChannels;

    % Read in the file
    reader.setSeries(tileIdx - 1);
    planes = reader.getImageCount;
    StackResults(tileIdx).numZ = planes/numChannels;
    Atest = bfGetPlane(reader, 1);
    ASz = size(Atest);
%     szIncX = floor(ASz(2)/numSegments);
%     szIncY = floor(ASz(1)/numSegments);
    
%*************************************************************************%
%**                           FIND THE ROIs                             **%
%*************************************************************************%
    if restartAnalysis || ~isfield(StackResults(tileIdx), 'numROIs') || isempty(StackResults(tileIdx).numROIs) || StackResults(tileIdx).numROIs <= 0
        clear ASum ASumDisp
        % Prepare the summed intensity projection image
        if restartAnalysis || ~isfield(StackResults(tileIdx), 'sumPrepared') || isempty(StackResults(tileIdx).sumPrepared)
            if StackResults(tileIdx).segChannel
                chList = StackResults(tileIdx).segChannel;
            else
                chList = 1:numChannels;
            end
            chCount = 0;
            for chIdx = chList
                chCount = chCount+1;
                ASum(:,:, chCount) = double(zeros(ASz)); 
                    % Loop through Z
                    for zIdx = 1:StackResults(tileIdx).numZ
                        disp(['Ch. ', num2str(chIdx), ', Z ', num2str(zIdx)])
                        % Read in the image
                        iPlane = reader.getIndex(zIdx-1, chIdx-1, 0) + 1;
                        A = bfGetPlane(reader, iPlane);
                        ASum(:,:, chCount) = ASum(:,:, chCount) + double(A);
                    end
                ASumDisp(:,:, chCount) = sgautoscale(ASum(:,:, chCount));
            end

            if size(ASumDisp, 3) < 2
                ASumDisp(:,:, 2) = double(zeros(ASz)); 
                ASum(:,:, 2) = double(zeros(ASz)); 
            end
            if size(ASumDisp, 3) < 3
                ASumDisp(:,:, 3) = double(zeros(ASz)); 
                ASum(:,:, 3) = double(zeros(ASz)); 
            end
            StackResults(tileIdx).sumPrepared = 1;
            StackResults(tileIdx).ASum = mean(double(ASum),3);
            StackResults(tileIdx).ASumDisp = mean(double(ASumDisp),3);
            save(dataOutputFilename, 'StackResults', 'ASumDisp');
        end

        h = figure(300); 
        pos = get(h, 'Position'); 
        imshow((ASumDisp(:,:,1:3)));
        set(h, 'Position', pos);
        
        % automatically find the ROIs
        % initialize the fields that segmentImage will assign
        StackResults(tileIdx).centroids = 0;
        StackResults(tileIdx).ROILocX = 0;
        StackResults(tileIdx).ROILocY = 0;
        StackResults(tileIdx).numROIs = 0;
            
        tmp = segmentImage(StackResults(tileIdx));
        StackResults(tileIdx) = tmp;
        numROIs = StackResults(tileIdx).numROIs;
        
        % delete ROIs that are touching the edge of the image
        deleteROI = zeros(numROIs,1);
        for roiIdx = 1:numROIs
            % find the distance from the ROI edge to the closest edge of
            % the image
            dist2edge = min([StackResults(tileIdx).ROILocX(roiIdx), (ASz(1)-(StackResults(tileIdx).ROILocX(roiIdx) + roiSz)), StackResults(tileIdx).ROILocY(roiIdx), (ASz(2)-(StackResults(tileIdx).ROILocY(roiIdx) + roiSz))]);
            % if the distance is smaller than 1 (no pixels between the ROI and the image edge), flag the ROI for deletion 
            if dist2edge < 1
                deleteROI(roiIdx) = 1;
            end
        end
        % delete the flagged ROIs
        StackResults(tileIdx).ROILocX(deleteROI ~=0) = [];
        StackResults(tileIdx).ROILocY(deleteROI ~=0) = [];
        if length(StackResults(tileIdx).ROILocX) < 1
            StackResults(tileIdx).ROILocX = [];
            StackResults(tileIdx).ROILocY = [];
        end
        StackResults(tileIdx).centroids{1,1}((deleteROI ~=0),:) = [];
        StackResults(tileIdx).numROIs = StackResults(tileIdx).numROIs-sum(deleteROI);
        disp(['Removed ', num2str(sum(deleteROI)), ' ROIs touching the frame edge']);
  
        if StackResults(tileIdx).numROIs
            % store the 'absolute ROI positions' -- including the stage
            % position data.  ROIabs = ROI*XYvoxelSz + stagePosition
            StackResults(tileIdx).AbsROILocX = StackResults(tileIdx).ROILocX.*StackResults(tileIdx).voxelSizeXYdefaultValue.doubleValue + StackResults(tileIdx).stageX;
            StackResults(tileIdx).AbsROILocY = StackResults(tileIdx).ROILocY.*StackResults(tileIdx).voxelSizeXYdefaultValue.doubleValue + StackResults(tileIdx).stageY;

            % store the tile number and ROI index for each ROI
            StackResults(tileIdx).ROITile = tileIdx.*ones(length(StackResults(tileIdx).ROILocX),1);
            StackResults(tileIdx).ROIIdx = (1:length(StackResults(tileIdx).ROILocX))';
        end
        save(dataOutputFilename, 'StackResults', 'ASumDisp');
    end
    
end
% now delete the summed images to save space and then save StackResults
[StackResults(1:numTiles).ASum] = deal([]);
[StackResults(1:numTiles).ASumDisp] = deal([]);
save(dataOutputFilename, 'StackResults', 'ASumDisp');

%% 


%*************************************************************************%
%**       FILTER THE ROIs TO FIND DUPLICATES (IN OVERLAP REGION)        **%
%*************************************************************************%
% initialize the CompiledResults ROIs table, which stores information 
% on all of the ROIs.  Columns of this table include the (1) the absolute
% X position, (2) the absolute Y position, (3) the in-image X position, (4)
% the in-image Y position, (5) the tile number in which the ROI is
% situated, (6) the ROI index in the original tile data structure for that
% ROI
X = cat(1,StackResults.ROILocX);
Y = cat(1,StackResults.ROILocY);
absX = cat(1,StackResults.AbsROILocX);
absY = cat(1,StackResults.AbsROILocY);
tile = cat(1,StackResults.ROITile);
inTileROIIdx = cat(1,StackResults.ROIIdx);
CompiledResults = table(X,Y,absX,absY,tile,inTileROIIdx);
figure(300); clf; plot(absX, absY, '.'); title('Absolute XY positions of ROIs before duplicate removal');



% now that all of the data are compiled, find table entries with the same
% absolute position within a tolerance for X and Y
% the tolerance for saying that two ROIs are at the same location is whether 
% the difference in their X and Y coordinates is smaller than the absolute 
% ROI size (microns)
tolerance = roiSz*StackResults(tileIdx).voxelSizeXYdefaultValue.doubleValue; % set the tolerance
[C,IA,IC] = uniquetol([absX,absY], tolerance, 'DataScale', 1, 'ByRows', true); % find the unique rows within the tolerance
% get a list of all the duplicate ROIs
uniqueIndices = unique(IC);
numRepetitions = histc(IC, uniqueIndices);
repeatingIndices = uniqueIndices(numRepetitions>1);

% initialize a variable to keep track of all of the tiles we need to go
% back and delete ROIs from 
tilesToModify = [];

% loop through all of the sets of duplicates
for dupSetIdx = 1:length(repeatingIndices)
    % find the duplicate entries for this set
    dupEntries = find(IC == repeatingIndices(dupSetIdx));
    disp('Found duplicate ROIs:')
    
%     % throw an error if the same ROI is acquired more than 4 times (should not
%     % be possible when tiling -- an ROI should only be able to overlap 4 fields of view)
%     if length(dupEntries) > 4
%         error(['Found ', num2str(length(dupEntries)), ' duplicates of the same ROI.  Check the code - this should not be possible with tiled acquisition.'])
%     end
    disp(CompiledResults(dupEntries, :));
    
    % For each duplicate, calculate the distance from the edge of the image
    % to the ROI centre
    clear dist2edge;
    for dupIdx=1:length(dupEntries)
        dupVal = dupEntries(dupIdx);
        dist2edge(dupIdx) = min([X(dupVal),(ASz(1)-(X(dupVal) + roiSz)), ...
            Y(dupVal), (ASz(2)-(Y(dupVal) + roiSz))]);
    end
    [tmp, bestIdx] = max(dist2edge);
    roiToKeep = dupEntries(bestIdx);
    disp('Keeping:')
    disp(CompiledResults(roiToKeep, :));
    
    % throw out all duplicates except the one in which the ROI is furthest
    % from the image edge
    for dupIdx=1:length(dupEntries)
        dupVal = dupEntries(dupIdx);
        if dupVal ~= roiToKeep
            % plot this point in red to show we have removed it
            figure(300); hold on; plot(StackResults(tile(dupVal)).AbsROILocX(inTileROIIdx(dupVal)), StackResults(tile(dupVal)).AbsROILocY(inTileROIIdx(dupVal)), '.r');
            
            % flag the duplicate ROIs for removal from the individual tile 
            % structure members
            tilesToModify = unique([tilesToModify,tile(dupVal)]);
            StackResults(tile(dupVal)).ROILocX(inTileROIIdx(dupVal)) = DELETEFLAG;
            StackResults(tile(dupVal)).ROILocY(inTileROIIdx(dupVal)) = DELETEFLAG;
            StackResults(tile(dupVal)).AbsROILocX(inTileROIIdx(dupVal)) = DELETEFLAG;
            StackResults(tile(dupVal)).AbsROILocY(inTileROIIdx(dupVal)) = DELETEFLAG;
            StackResults(tile(dupVal)).ROITile(inTileROIIdx(dupVal)) = DELETEFLAG;
            StackResults(tile(dupVal)).numROIs = StackResults(tile(dupVal)).numROIs-1;
        end
    end
end

% remove the duplicates from the individual tile structure members
for tileIdx = 1:length(tilesToModify)
    tileVal = tilesToModify(tileIdx);
    StackResults(tileVal).ROILocX(StackResults(tileVal).ROILocX==DELETEFLAG) = [];
    StackResults(tileVal).ROILocY(StackResults(tileVal).ROILocY==DELETEFLAG) = [];
    StackResults(tileVal).AbsROILocX(StackResults(tileVal).AbsROILocX==DELETEFLAG) = [];
    StackResults(tileVal).AbsROILocY(StackResults(tileVal).AbsROILocY==DELETEFLAG) = [];
    StackResults(tileVal).ROITile(StackResults(tileVal).ROITile==DELETEFLAG) = [];
    StackResults(tileVal).ROIIdx = (1:length(StackResults(tileVal).ROILocX))';
    
    if ~isequal(length(StackResults(tileVal).ROILocX), length(StackResults(tileVal).ROILocY), length(StackResults(tileVal).AbsROILocX), ...
            length(StackResults(tileVal).AbsROILocY), length(StackResults(tileVal).ROITile), length(StackResults(tileVal).ROIIdx), ...
            StackResults(tileVal).numROIs)
        warning('Deletion of duplicate ROIs yielded different numbers of ROIs for different structure elements.');
        disp(['ROILocX: ', num2str(length(StackResults(tileVal).ROILocX)), ...
            ', ROILocY: ', num2str(length(StackResults(tileVal).ROILocY)), ...
            ', AbsROILocX: ', num2str(length(StackResults(tileVal).AbsROILocX)), ...
            ', AbsROILocY: ', num2str(length(StackResults(tileVal).ROILocY)), ...
            ', ROITile: ', num2str(length(StackResults(tileVal).ROITile)), ...
            ', ROIIdx: ', num2str(length(StackResults(tileVal).ROIIdx)), ...
            ', numROIs: ', num2str(StackResults(tileVal).numROIs)]);
    end
end

% remake the compiled table of ROIs
X = cat(1,StackResults.ROILocX);
Y = cat(1,StackResults.ROILocY);
absX = cat(1,StackResults.AbsROILocX);
absY = cat(1,StackResults.AbsROILocY);
tile = cat(1,StackResults.ROITile);
inTileROIIdx = cat(1,StackResults.ROIIdx);
CompiledResults = table(X,Y,absX,absY,tile,inTileROIIdx);

h = figure(300);
saveas(h, [StackResults(tileIdx).outputFilenameBase, 'RoiLocsPreDeletion.fig']);

figure(301); clf; plot(absX, absY, '.'); title('Absolute XY positions of ROIs after duplicate removal');
saveas(h, [StackResults(tileIdx).outputFilenameBase, 'RoiLocsPostDeletion.fig']);

save(dataOutputFilename, 'StackResults', 'ASumDisp', 'CompiledResults');

%% 


                
%*************************************************************************%
%**          LOOP THROUGH THE DATA TO SUM INTENSITY PROFILES            **%
%*************************************************************************%
% Loop through the list of stacks
for tileIdx = 1:numTiles
%for tileIdx = 11:13
    if ~isfield(StackResults(tileIdx), 'IsumBkgSub') || isempty(StackResults(tileIdx).IsumBkgSub) || ~any(StackResults(tileIdx).IsumBkgSub(:) > 0)
        reader.setSeries(tileIdx - 1);
        planes = reader.getImageCount;
        
        
        disp(['Tile ', num2str(tileIdx)])
        % Set up background region for each ROI
        for roiIdx = 1:StackResults(tileIdx).numROIs
            [xStart, yStart, xEnd, yEnd] = findROILims(StackResults(tileIdx), roiIdx, roiSz, ASz, 1);

            % choose background region as the pixels surrounding
            % each ROI (set by wBkg), with small region to calculate
            % background-subtracted noise in between
            roiRegion = logical(zeros(ASz));
            roiRegion(xStart:xEnd, yStart:yEnd) = 1;
            roiRegionDilNoise = imdilate(roiRegion, strel('disk', wNoiseBkg));
            roiNoise(:,:, roiIdx) = roiRegionDilNoise-roiRegion;
            roiRegionDilBkg = imdilate(roiRegionDilNoise, strel('disk', wBkg));
            roiBkg(:,:, roiIdx) = roiRegionDilBkg-roiRegionDilNoise;
        end
        
        % Find startIdx for each ROI
        StackResults(tileIdx).startIdx = zeros(StackResults(tileIdx).numROIs,1);
        done = 0;
        zIdx = 0;
        if StackResults(tileIdx).numROIs
            % Set up the list of ROIs
            roiList = 1:StackResults(tileIdx).numROIs;
            while ~done % looping, increasing zNum each loop
                zIdx = zIdx + 1;
                disp(['Finding topZs: Z=', num2str(zIdx)])
                ASum = double(zeros(ASz));
                % create summed image of all colour channels
                for chIdx = chList
                    iPlane = reader.getIndex(StackResults(tileIdx).numZ-(zIdx), chIdx-1, 0) + 1;
                    A = bfGetPlane(reader, iPlane);
                    ASum = ASum + double(sgautoscale(A));
                end
                % normalize and segment the image
                ASumNorm = double(ASum)./double(max(max(ASum)));
                ABW = imbinarize(ASumNorm, 'adaptive', 'sensitivity', 0.9*fudge);
                se = strel('disk', 5);
                ABW = imclose(ABW, se);
                se = strel('disk', 7);
                ABW = imopen(ABW, se);
                ABW = imfill(ABW, 'holes');
                ABW = imclose(ABW, se);
                se = strel('disk', 10);
                ABW = imopen(ABW, se);
                ABW = imclearborder(ABW);
                ABW = bwareaopen(ABW, minSpotSize);
                figure(250); clf; imshow(sgautoscale(ASumNorm)); hold on; 
                plot(StackResults(tileIdx).ROILocX+roiSz/2, StackResults(tileIdx).ROILocY+roiSz/2, 'xr'); hold on;
                plot(StackResults(tileIdx).ROILocX(StackResults(tileIdx).startIdx ~=0)+roiSz/2, StackResults(tileIdx).ROILocY(StackResults(tileIdx).startIdx ~=0)+roiSz/2, 'xg'); pause(0.01);
                % loop through the ROIs to see if we see protein spots of the
                % appropriate size
                for roiIdx = roiList
                    [xStart, yStart, xEnd, yEnd] = findROILims(StackResults(tileIdx), roiIdx, roiSz, ASz, 1);
                    ABWROI = ABW(xStart:xEnd, yStart:yEnd);
                    % if we see the appropriate size spot, store this ZIdx as
                    % startIdx for this ROI, and remove it from the list of
                    % ROIs so we don't loop through it anymore
                    if(nnz(ABWROI)) > minSpotSize
                        disp(['ROI ', num2str(roiIdx), ' startIdx = ', num2str(zIdx)]); 
                        StackResults(tileIdx).startIdx(roiIdx) = zIdx;
                        % remove this ROI from the list
                        roiList = roiList(roiList ~= roiIdx);
                    end
                end
                % have we found all of the start indices (for all ROIs)?
                % If not, keep looping
                if all(StackResults(tileIdx).startIdx) || zIdx >= StackResults(tileIdx).numZ/2
                    done = 1;
                end
            end



            % Loop through the colour channels
            for chIdx = 1:numChannels
                    % Loop through Z
                    for zIdx = 1:StackResults(tileIdx).numZ
                        disp(['Ch. ', num2str(chIdx), ', Z ', num2str(zIdx)])
                        StackResults(tileIdx).Z(zIdx) = double(StackResults(tileIdx).voxelSizeZdefaultValue)*(zIdx-1);
                        % Read in the image
                        iPlane = reader.getIndex((zIdx-1), chIdx-1, 0) + 1;
                        A = bfGetPlane(reader, iPlane);
                        ASz = size(A);

                        % Loop through the ROIs
                        for roiIdx = 1:StackResults(tileIdx).numROIs
                            [xStart, yStart, xEnd, yEnd] = findROILims(StackResults(tileIdx), roiIdx, roiSz, ASz, 1);


                            % get the average background value for this ROI
                            StackResults(tileIdx).avgBkg(chIdx, roiIdx, zIdx) = mean2(A(roiBkg(:,:, roiIdx)~=0));
                            
                            % get the average value in noise calculation
                            % region for this ROI
                            StackResults(tileIdx).avgBkgNoiseRegion(chIdx, roiIdx, zIdx) = mean2(A(roiNoise(:,:, roiIdx)~=0));
                            StackResults(tileIdx).avgBkgNoiseRegionBkgSub(chIdx, roiIdx, zIdx) = mean2(A(roiNoise(:,:, roiIdx)~=0)) ...
                                - StackResults(tileIdx).avgBkg(chIdx, roiIdx, zIdx);

                            % take the region
                            Aregion = A(xStart:xEnd, yStart:yEnd);

                            % Sum all the intensities in the ROI and store
                            StackResults(tileIdx).Isum(chIdx, roiIdx, zIdx) = sum(double(Aregion(:)));
                            StackResults(tileIdx).IsumBkgSub(chIdx, roiIdx, zIdx) = ...
                                sum(double(Aregion(:)))-length(Aregion(:))*double(StackResults(tileIdx).avgBkg(chIdx, roiIdx, zIdx));

                        end % end ROI loop
                    end % End Z loop
                    % flip the Z stack if it was acquired from bottom to top of the gel
                    if StackResults(tileIdx).flipStack
                        StackResults(tileIdx).Z = flip(StackResults(tileIdx).Z);
                    end

                    % compute the noise as the standard deviation of the
                    % background intensity with Z)
                    StackResults(tileIdx).bkgNoiseMagnitude(chIdx, roiIdx) = std(StackResults(tileIdx).avgBkg(chIdx, roiIdx, :));

            end % End colour loop
        end

        save(dataOutputFilename, 'StackResults', 'ASumDisp', 'CompiledResults');
    end % End stack loop
end % end check loop

for tileIdx = 1:numTiles
%for tileIdx = 24:26
    if StackResults(tileIdx).numROIs
         h = figure(tileIdx); clf;
        for roiIdx = 1:StackResults(tileIdx).numROIs
            % use the TopZ (startIdx) value as the relative Z location
            StackResults(tileIdx).ZScaled(roiIdx, :) = StackResults(tileIdx).Z-(double(StackResults(tileIdx).startIdx(roiIdx))*double(StackResults(tileIdx).voxelSizeZdefaultValue));
            for chIdx = 1:numChannels
                % Calculate the normalized and minimum-subtracted profile and plot 
                % StackResults(tileIdx).IsumMinSub(chIdx, roiIdx, :) = ...
                %    (StackResults(tileIdx).Isum(chIdx, roiIdx, :)-min(StackResults(tileIdx).Isum(chIdx, roiIdx, :)));
                h = figure(tileIdx); hold on;
                plot(squeeze(StackResults(tileIdx).ZScaled(roiIdx, :)), squeeze(StackResults(tileIdx).IsumBkgSub(chIdx, roiIdx, :)).*channelScaling(chIdx), 'Color', (defaultColours(chIdx,:)), 'Linewidth', 1);  
                xlim([min(StackResults(tileIdx).ZScaled(roiIdx, :)), max(StackResults(tileIdx).ZScaled(roiIdx, :))])
            end
        end

        % Add legend, x, and y labels
        xlabel('Z depth (\mum)')
        ylabel('Normalized fluor. intensity (AFU)')
        set(h, 'Position', [150 600 250 140]);

        % Set title
        title(StackResults(tileIdx).title);

        % Save
        saveas(h, [StackResults(tileIdx).outputFilenameBase, 'tile', num2str(tileIdx), '_IProfilesNoLeg.fig']);
        saveas(h, [StackResults(tileIdx).outputFilenameBase, 'tile', num2str(tileIdx), '_IProfilesNoLeg.svg']);
        legend(proteins, 'Location', 'Northeast')
        saveas(h, [StackResults(tileIdx).outputFilenameBase, 'tile', num2str(tileIdx), 'IProfiles.fig']);
        saveas(h, [StackResults(tileIdx).outputFilenameBase, 'tile', num2str(tileIdx), 'IProfiles.svg']);
    end

end

%% 

% Gaussian fit the peaks and extract the migration distances and peak
% widths
for tileIdx = 1:numTiles
    if ~isfield(StackResults(tileIdx), 'intensityFitBounds') || isempty(StackResults(tileIdx).intensityFitBounds) || ~any(any(StackResults(tileIdx).intensityFitBounds > 0))
    % Choose the start of the gel:
        if StackResults(tileIdx).numROIs
            h = figure(11); clf;
            for roiIdx = 1:StackResults(tileIdx).numROIs
                % use the TopZ (startIdx) value as the relative Z location
                StackResults(tileIdx).ZScaled(roiIdx, :) = StackResults(tileIdx).Z-(double(StackResults(tileIdx).startIdx(roiIdx))*double(StackResults(tileIdx).voxelSizeZdefaultValue));
                for chIdx = 1:numChannels
                    % Calculate the normalized and minimum-subtracted profile and plot 
                    % StackResults(tileIdx).IsumMinSub(chIdx, roiIdx, :) = ...
                    %    (StackResults(tileIdx).Isum(chIdx, roiIdx, :)-min(StackResults(tileIdx).Isum(chIdx, roiIdx, :)));
                    h = figure(11); hold on;
                    plot(squeeze(StackResults(tileIdx).Z), squeeze(StackResults(tileIdx).IsumBkgSub(chIdx, roiIdx, :)).*channelScaling(chIdx), 'Color', (defaultColours(chIdx,:)), 'Linewidth', 1);  
                    xlim([min(StackResults(tileIdx).Z), max(StackResults(tileIdx).ZScaled(roiIdx, :))])
                end
            end
            title(['Stack ', num2str(tileIdx), ': click on the start of the gel or press n to disregard all ROIs from this tile']);
            [x,y,button] = ginput(1);

            % remove this tile from analysis if none of its ROIs seem usable
            if button == 'n' || button == 'N'
                StackResults(tileIdx).numROIs = 0;
            end

            for roiIdx = 1:StackResults(tileIdx).numROIs
                StackResults(tileIdx).startIdx(roiIdx) = uint16(x/double(StackResults(tileIdx).voxelSizeZdefaultValue));
                StackResults(tileIdx).ZScaled(roiIdx, :) = StackResults(tileIdx).Z-x;
            end
        end
    end
end
    
%% 

for tileIdx = 1:numTiles    
        %for tileIdx = 24:26
    if StackResults(tileIdx).numROIs
        if ~isfield(StackResults(tileIdx), 'intensityFitBounds') || isempty(StackResults(tileIdx).intensityFitBounds) || ~any(any(StackResults(tileIdx).intensityFitBounds > 0))

            % fit the peaks and measure the migration distances and peak widths 
                % for the protein bands
                % set up some of the fields that the measuremigration function needs
                StackResults(tileIdx).stackIdx = tileIdx;
                StackResults(tileIdx).peaknames = proteins;

                % initialize the fields that measuremigration will assign
                StackResults(tileIdx).intPkFitA = 0;
                StackResults(tileIdx).intPkFitMu = 0;
                StackResults(tileIdx).intPkFitSigma = 0;
                StackResults(tileIdx).intPkFitRSq = 0;
                StackResults(tileIdx).intensityFitBounds = 0;
                StackResults(tileIdx).intPkSNR = 0;
                StackResults(tileIdx).intPkSNRReg = 0;

                %StackResults(tileIdx) = measuremigration(StackResults(tileIdx));
                tmp = measuremigrationLS(StackResults(tileIdx));
                StackResults(tileIdx) = tmp;
                save(dataOutputFilename, 'StackResults', 'ASumDisp', 'CompiledResults');
        end
    end
end


%% 


% % collect all of the migration distance data for the same experimental
% % conditions
%close(reader);
% % find all of the unique experimental conditions from the title field
% uniqueConditions = unique({StackResults(~isempty(StackResults.title)).title});
% numConditions = length(uniqueConditions);
numConditions = 1;
% % loop through the unique experimental conditions, collecting the data for
% % each stack of that condition
for condIdx = 1:numConditions
    % find all the stacks of this experimental condition
    matchingStacks = find(strcmp({StackResults(11).title}, {StackResults.title})>0);
    
    % loop through the matching stacks and collect the data
    for chIdx = 1:numChannels
        for pkIdx = 1:numPeaks
            collectedMu = [];
            collectedSigma = [];
            collectedXY = [];
            for tileIdx = matchingStacks
                if StackResults(tileIdx).numROIs
                    collectedMu = [collectedMu;(nonzeros(squeeze(StackResults(tileIdx).intPkFitMu(chIdx,pkIdx,:)))-198)];
                    collectedSigma = [collectedSigma;nonzeros(squeeze(StackResults(tileIdx).intPkFitSigma(chIdx,pkIdx,:)))];
                    %collectedXY = [collectedXY;nonzeros(squeeze(StackResults(tileIdx).intPkFitXYRes(chIdx,pkIdx,:)))];
                end
            end
            migrationCollected(condIdx, chIdx) = {collectedMu};
            sigmaCollected(condIdx, chIdx) = {collectedSigma};
            xyResCollected(condIdx, chIdx) = {collectedXY};
            
            migrationMean(condIdx, chIdx) = mean(collectedMu);
            migrationErr(condIdx, chIdx) = std(collectedMu);
            
            sigmaMean(condIdx, chIdx) = mean(collectedSigma);
            collSigmaSq = collectedSigma.^2;
            sigmaSqMean(condIdx, chIdx) = mean(collSigmaSq);
            sigmaErr(condIdx, chIdx) = std(collectedSigma);
            sigmaSqErr(condIdx, chIdx) = std(collSigmaSq);
            
%             xyResMean(condIdx, chIdx) = mean(collectedXY);
%             xyResSq = collectedXY.^2;
%             xyResSqMean(condIdx, chIdx) = mean(xyResSq);
%             xyResErr(condIdx, chIdx) = std(collectedXY);
%             xyResSqErr(condIdx, chIdx) = std(xyResSq);
            
        end
    end
    

end

h = figure; plotSpread(migrationCollected, 'distributionColors', [{'m'},{'c'}], 'xNames', [{'GAPDH'},{'Actinin'}]);
ylabel('Migration distance (\mum)');
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_migrationBeeswarm.fig']);
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_migrationBeeswarm.svg']);

h = figure; plotSpread(sigmaCollected, 'distributionColors', [{'m'},{'c'}], 'xNames', [{'GAPDH'},{'Actinin'}]);
ylabel('Peak width (\mum)');
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_pkWBeeswarmCorrected.fig']);
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_pkWBeeswarmCorrected.svg']);

%CompiledResults.Protein1Migration = cat(squeeze(StackResults.intPkFitMu(1,1,:)));
%CompiledResults.Protein2Migration = cat(squeeze(StackResults.intPkFitMu(2,1,:)));
save(dataOutputFilename, 'StackResults', 'migrationCollected', 'sigmaCollected', 'xyResCollected', 'CompiledResults', 'ASumDisp');


%% 
for condIdx = 1:numConditions
    % find all the stacks of this experimental condition
    matchingStacks = find(strcmp({StackResults(11).title}, {StackResults.title})>0);
    
    for tileIdx = matchingStacks
        for roiIdx = 1:numROIs
                if StackResults(tileIdx).numROIs
                    collectedMu = [collectedMu;(nonzeros(squeeze(StackResults(tileIdx).intPkFitMu(chIdx,pkIdx,:)))-186)];
                    collectedSigma = [collectedSigma;nonzeros(squeeze(StackResults(tileIdx).intPkFitSigma(chIdx,pkIdx,:)))];
                    %collectedXY = [collectedXY;nonzeros(squeeze(StackResults(tileIdx).intPkFitXYRes(chIdx,pkIdx,:)))];
                end
        end
    % loop through the matching stacks and collect the data
    for chIdx = 1:numChannels
        for pkIdx = 1:numPeaks
            collectedMu = [];
            collectedSigma = [];
            collectedXY = [];
            
            migrationCollected(condIdx, chIdx) = {collectedMu};
            sigmaCollected(condIdx, chIdx) = {collectedSigma};
            xyResCollected(condIdx, chIdx) = {collectedXY};
            
            migrationMean(condIdx, chIdx) = mean(collectedMu);
            migrationErr(condIdx, chIdx) = std(collectedMu);
            
            sigmaMean(condIdx, chIdx) = mean(collectedSigma);
            collSigmaSq = collectedSigma.^2;
            sigmaSqMean(condIdx, chIdx) = mean(collSigmaSq);
            sigmaErr(condIdx, chIdx) = std(collectedSigma);
            sigmaSqErr(condIdx, chIdx) = std(collSigmaSq);
            
%             xyResMean(condIdx, chIdx) = mean(collectedXY);
%             xyResSq = collectedXY.^2;
%             xyResSqMean(condIdx, chIdx) = mean(xyResSq);
%             xyResErr(condIdx, chIdx) = std(collectedXY);
%             xyResSqErr(condIdx, chIdx) = std(xyResSq);
        end
            
        end
    end
    

end

h = figure; plotSpread(migrationCollected, 'distributionColors', [{'m'},{'c'}], 'xNames', [{'GAPDH'},{'Actinin'}]);
ylabel('Migration distance (\mum)');
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_migrationBeeswarmCorrected.fig']);
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_migrationBeeswarmCorrected.svg']);

h = figure; plotSpread(sigmaCollected, 'distributionColors', [{'m'},{'c'}], 'xNames', [{'GAPDH'},{'Actinin'}]);
ylabel('Peak width (\mum)');
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_pkWBeeswarmCorrected.fig']);
saveas(h, [StackResults(tileIdx).outputFilenameBase, '_pkWBeeswarmCorrected.svg']);

%CompiledResults.Protein1Migration = cat(squeeze(StackResults.intPkFitMu(1,1,:)));
%CompiledResults.Protein2Migration = cat(squeeze(StackResults.intPkFitMu(2,1,:)));
save(dataOutputFilename, 'StackResults', 'migrationCollected', 'sigmaCollected', 'xyResCollected', 'CompiledResults', 'ASumDisp');





% tEP = [10,20];
% tDiff = [20,30];
% hMu = figure; hold on;
% hSigma = figure; hold on;
% hXY = figure; hold on;
% for chIdx = 1:numChannels
%     muPlot = [];
%     sigmaPlot = [];
%     xyPlot = [];
%     muX = [];
%     sigmaX = [];
%     xyX = [];
%     for condIdx = 1:size(migrationCollected, 1)
%         muPlot = [muPlot;migrationCollected{condIdx, chIdx}];
%         sigmaPlot = [sigmaPlot;(sigmaCollected{condIdx, chIdx})];
%         xyPlot = [xyPlot;(xyResCollected{condIdx, chIdx})];
% %         sigmaPlot = [sigmaPlot;(sigmaCollected{condIdx, chIdx}).^2];
% %         xyPlot = [xyPlot;(xyResCollected{condIdx, chIdx}).^2];
%         
%         muX = [muX; tEP(condIdx).*ones(length(migrationCollected{condIdx, chIdx}), 1)];
%         sigmaX = [sigmaX; tDiff(condIdx).*ones(length(sigmaCollected{condIdx, chIdx}), 1)];
%         xyX = [xyX; tDiff(condIdx).*ones(length(xyResCollected{condIdx, chIdx}), 1)];
%     end
%         
%     figure(hMu); plot(muX, muPlot, ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('EP time (s)'); ylabel('Migration distance (\mum)');
%     xlim([0, 25]);
%     figure(hSigma); plot(sigmaX, sigmaPlot, ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('Diffusion time (s)'); ylabel('\sigma (\mum)');
%     %xlabel('Diffusion time (s)'); ylabel('\sigma^2 (\mum^2)');
%     xlim([0, 35]);
%     figure(hXY); plot(xyX, xyPlot, ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('Diffusion time (s)'); ylabel('XY peak width (\mum)');
%     %xlabel('Diffusion time (s)'); ylabel('Squared XY peak width (\mum^2)');
%     xlim([0, 35]);
% 
% %     figure(hMu); errorbar(tEP, migrationMean(:, chIdx), migrationErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('EP time (s)'); ylabel('Migration distance (\mum)');
% %     xlim([0, 25]);
% %     figure(hSigma); errorbar(tDiff, sigmaSqMean(:, chIdx), sigmaSqErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('Diffusion time (s)'); ylabel('\sigma^2 (\mum^2)');
% %     xlim([0, 35]);
% %     figure(hXY); errorbar(tDiff, xyResSqMean(:, chIdx), xyResSqErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx,:))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('EP time (s)'); ylabel('Squared XY peak width (\mum^2)');
% %     xlim([0, 35]);
% end
% figure(hMu); hold on;
% x = 0:25;
% plot(x, 30*x+1.5, '-.c');
% plot(x, 49*x-210, '-.m');
% plot(x, 24*x-75, '-.r');
% 
