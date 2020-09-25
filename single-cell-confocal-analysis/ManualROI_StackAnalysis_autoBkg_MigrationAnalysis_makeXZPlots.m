%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%**            ZEP STACK QUICK Z PROFILE ANALYSIS SCRIPT                **%
%**                         July 22ND, 2019                             **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

% This wrapper script does image processing analysis on a series of Z stacks 
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
else % MacOS path
    dataAnalysisRoot = 'path2';
end
addpath(dataAnalysisRoot);
if ispc
    addpath([dataAnalysisRoot, '\plotSpread']);
else
    addpath([dataAnalysisRoot, '/plotSpread']);
end

% This script requires the Open Microscopy Bio-Formats MATLAB toolbox,
% available: http://www.openmicroscopy.org/bio-formats/downloads/

% add the bioformats image import code to the path
if ispc
    bioformatsPath = 'path3';
else
    bioformatsPath = 'path4';
end
addpath(bioformatsPath);


%*************************************************************************%
%**                     SET UP THE CONFIGURATIONS                       **%
%*************************************************************************%

% Read in channel configuration information
% load in the configuration .csv file
configFilename = 'Initialization.xlsx';
dataOutputFilename = '2019-11-21_Analysis_BkgSub_DATA_MigrationAnalysis.mat';
configData = readtable(configFilename);
numStacks = size(configData, 1);
%numSegments = 4; % the number by which to divide the image in x and y to create ROIs (numSegments^2 ROIs)
roiSz = 300; % size in pixels of the length and width of the rectangular ROIs
wBkg = 20; % width of the background 'gutter' region surrounding each ROI
wNoiseBkg = 10; % width of the extra background 'gutter' region surrounding each ROI used for calculation of noise in background-subtracted intensity profiles
wBkgXY = 5; % portion of the X-Y ROI to use as background when measuring x and y intensity profiles for XY resolution calcs
maxNumROIs = 4; % maximum number of ROIs to analyze in each image
makeXZPlots = 1; % generate and save X-Z summed 2D images for each ROI
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
defaultColours = channelData{:, 2};

%*************************************************************************%
%**                     SET UP THE CONFIGURATIONS                       **%
%*************************************************************************%
% Loop through the list of stacks
for stackIdx = 1:numStacks
    disp(['Stack ', num2str(stackIdx)])
    % Set up the filenames
    StackResults(stackIdx).proteins = proteins;
    StackResults(stackIdx).defaultColours = defaultColours;
    StackResults(stackIdx).inputFilename = cell2mat(configData{stackIdx, 1});
    StackResults(stackIdx).outputFilenameBase = cell2mat(configData{stackIdx, 2}); 
    StackResults(stackIdx).title = cell2mat(configData{stackIdx, 3});
    StackResults(stackIdx).flipStack = configData{stackIdx, 4};
    StackResults(stackIdx).topZ = configData{stackIdx, 5};
    StackResults(stackIdx).roiSz = roiSz;
    StackResults(stackIdx).wBkgXY = wBkgXY;
    
    numPeaks = size(StackResults(stackIdx).proteins, 2);

    % Get the XY and Z resolution
    StackResults(stackIdx).imgReader = bfGetReader(StackResults(stackIdx).inputFilename); % use the bioformats toolbox to create the image reader for this stack
    reader = StackResults(stackIdx).imgReader;
    data = reader.getMetadataStore();
    StackResults(stackIdx).voxelSizeXYdefaultValue = data.getPixelsPhysicalSizeX(0).value();
    StackResults(stackIdx).voxelSizeZdefaultValue = data.getPixelsPhysicalSizeZ(0).value();

    % Check that we have the right number of colour channels
    if data.getChannelCount(0) ~= numChannels
        warning('The number of defined colour channels is not equal to the image colour channels');
        numChannels = min([data.getChannelCount(0), numChannels]);
    end
    StackResults(stackIdx).numChannels = numChannels;

    % Read in the file
    reader.setSeries(0);
    planes = reader.getImageCount;
    StackResults(stackIdx).numZ = planes/data.getChannelCount(0);
    Atest = bfGetPlane(reader, 1);
    ASz = size(Atest);
%     szIncX = floor(ASz(2)/numSegments);
%     szIncY = floor(ASz(1)/numSegments);
    
%*************************************************************************%
%**                 ALLOW THE USER TO SELECT THE ROIs                   **%
%*************************************************************************%
    if ~isfield(StackResults(stackIdx), 'numROIs') || isempty(StackResults(stackIdx).numROIs) || StackResults(stackIdx).numROIs <= 0
        % Prepare the summed intensity projection image
        for chIdx = 1:numChannels
            ASum(:,:, chIdx) = double(zeros(ASz)); 
                % Loop through Z
                for zIdx = 1:StackResults(stackIdx).numZ
                    disp(['Ch. ', num2str(chIdx), ', Z ', num2str(zIdx)])
                    % Read in the image
                    iPlane = reader.getIndex(zIdx-1, chIdx-1, 0) + 1;
                    A = bfGetPlane(reader, iPlane);
                    ASum(:,:, chIdx) = ASum(:,:, chIdx) + double(A);
                end
            ASumDisp(:,:, chIdx) = sgautoscale(ASum(:,:, chIdx));
        end

        if size(ASumDisp, 3) < 2
            ASumDisp(:,:, 2) = double(zeros(ASz)); 
        end
        if size(ASumDisp, 3) < 3
            ASumDisp(:,:, 3) = double(zeros(ASz)); 
        end

        h = figure(300); 
        pos = get(h, 'Position'); 
        imshow((ASumDisp(:,:,1:3)));
        set(h, 'Position', pos);
        done = 0;

        while ~done
        % Find the ROIs for the image.  Choose a width specified by
            % EPData.ROIWidth, and a height of the full image height.
            % add the first ROI
            title(['Image ', num2str(stackIdx), ': click in the centre of the ROI, move if necessary, and then double-click on it.']);
            [x,y,button] = ginput(1)
            hRect = imrect(gca, [x-roiSz/2, y-roiSz/2, roiSz, roiSz]);
            position = wait(hRect); % position has the form [xmin, ymin, width, height]
            StackResults(stackIdx).ROILocX(1) = position(1);
            StackResults(stackIdx).ROILocY(1) = position(2);

            % add up to 3 more ROIs
            for i = 2:maxNumROIs
                %disp('testInLoop');
                title(['Image ', num2str(stackIdx), ': click to add another ROI.  Press y to proceed without adding more ROIs.']);
                button = 0;
                while (button ~= 1) && (button ~= 'y')
                    [x,y,button] = ginput(1);
                end
                if (button == 1)
                    title(['Image ', num2str(stackIdx), ': move ROI if necessary, and then double-click on it.']);
                    hRect = imrect(gca, [x-roiSz/2, y-roiSz/2, roiSz, roiSz]);
                    position = wait(hRect); % position has the form [xmin, ymin, width, height]
                    StackResults(stackIdx).ROILocX(i) = position(1);
                    StackResults(stackIdx).ROILocY(i) = position(2);
                    % save the number of ROIs that we have saved
                    StackResults(stackIdx).numROIs = i;
                end   
                if (button == 'y')
                    % save the number of ROIs that we have saved
                    StackResults(stackIdx).numROIs = i-1;
                    break;
                end
            end
            done = 1;
            save(dataOutputFilename, 'StackResults');
            
        end
    end
    
end
%% 

                
%*************************************************************************%
%**          LOOP THROUGH THE DATA TO SUM INTENSITY PROFILES            **%
%*************************************************************************%
% Loop through the list of stacks
for stackIdx = 1:numStacks
    if ~isfield(StackResults(stackIdx), 'IsumBkgSub') || isempty(StackResults(stackIdx).IsumBkgSub) || ~any(StackResults(stackIdx).IsumBkgSub(:) > 0)
        StackResults(stackIdx).imgReader = bfGetReader(StackResults(stackIdx).inputFilename); % use the bioformats toolbox to create the image reader for this stack
        reader = StackResults(stackIdx).imgReader;
        reader.setSeries(0);
        planes = reader.getImageCount;
        disp(['Stack ', num2str(stackIdx)])
        % Set up background region for each ROI
        for roiIdx = 1:StackResults(stackIdx).numROIs
            yStart = uint16(StackResults(stackIdx).ROILocX(roiIdx));
            yEnd = uint16(yStart+roiSz);
            xStart = uint16(StackResults(stackIdx).ROILocY(roiIdx));
            xEnd = uint16(xStart+roiSz);
            % check the values
            if xStart < 1
                xStart = xStart + (1-xStart);
                xEnd = xEnd + (1-xStart);
            end
            if yStart < 1
                yStart = yStart + (1-yStart);
                yEnd = yEnd + (1-yStart);
            end
            if xEnd > ASz(1)
                xEnd = xEnd - (xEnd - ASz(1));
                xStart = xStart - (xEnd - ASz(1));
            end
            if yEnd > ASz(2)
                yEnd = yEnd - (yEnd - ASz(2));
                yStart = yStart - (yEnd - ASz(2));
            end

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

        % Loop through the colour channels
        for chIdx = 1:numChannels
                % Loop through Z
                for zIdx = 1:StackResults(stackIdx).numZ
                    disp(['Ch. ', num2str(chIdx), ', Z ', num2str(zIdx)])
                    StackResults(stackIdx).Z(zIdx) = double(StackResults(stackIdx).voxelSizeZdefaultValue)*(zIdx-1);
                    % Read in the image
                    iPlane = reader.getIndex(zIdx-1, chIdx-1, 0) + 1;
                    A = bfGetPlane(reader, iPlane);
                    ASz = size(A);

                    % Loop through the ROIs
                    for roiIdx = 1:StackResults(stackIdx).numROIs
                        yStart = uint16(StackResults(stackIdx).ROILocX(roiIdx));
                        yEnd = uint16(yStart+roiSz);
                        xStart = uint16(StackResults(stackIdx).ROILocY(roiIdx));
                        xEnd = uint16(xStart+roiSz);
                        % check the values and adjust ROIs to be within image 
                        if xStart < 1
                            xStart = xStart + (1-xStart);
                            xEnd = xEnd + (1-xStart);
                        end
                        if yStart < 1
                            yStart = yStart + (1-yStart);
                            yEnd = yEnd + (1-yStart);
                        end
                        if xEnd > ASz(1)
                            xEnd = xEnd - (xEnd - ASz(1));
                            xStart = xStart - (xEnd - ASz(1));
                        end
                        if yEnd > ASz(2)
                            yEnd = yEnd - (yEnd - ASz(2));
                            yStart = yStart - (yEnd - ASz(2));
                        end


                        % get the average background value for this ROI
                        StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx) = mean2(A(roiBkg(:,:, roiIdx)~=0));

                        % get the average value in noise calculation
                        % region for this ROI
                        StackResults(stackIdx).avgBkgNoiseRegion(chIdx, roiIdx, zIdx) = mean2(A(roiNoise(:,:, roiIdx)~=0));
                        StackResults(stackIdx).avgBkgNoiseRegionBkgSub(chIdx, roiIdx, zIdx) = mean2(A(roiNoise(:,:, roiIdx)~=0)) ...
                            - StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx);

                        % take the region
                        Aregion = A(xStart:xEnd, yStart:yEnd);

                        % Sum all the intensities in the ROI and store
                        StackResults(stackIdx).Isum(chIdx, roiIdx, zIdx) = sum(double(Aregion(:)));
                        StackResults(stackIdx).IsumBkgSub(chIdx, roiIdx, zIdx) = ...
                            sum(double(Aregion(:)))-length(Aregion(:))*double(StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx));
                        
                        % Make the X-Z summed images
                        if makeXZPlots
                            StackResults(stackIdx).xzSumImg(chIdx, roiIdx, :, zIdx) = sum(double(Aregion));
                            StackResults(stackIdx).xzSumImgBkgSub(chIdx, roiIdx, :, zIdx) = sum(double(Aregion))-size(Aregion,1)*double(StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx));
                        end

%                         % get the average background value for this ROI
%                         StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx) = mean2(A(roiBkg(:,:, roiIdx)~=0));

%                         % take the region
%                         Aregion = A(xStart:xEnd, yStart:yEnd);
% 
%                         % Sum all the intensities in the image and store
%                         StackResults(stackIdx).Isum(chIdx, roiIdx, zIdx) = sum(double(Aregion(:)));
%                         StackResults(stackIdx).IsumBkgSub(chIdx, roiIdx, zIdx) = ...
%                             sum(double(Aregion(:)))-length(Aregion(:))*double(StackResults(stackIdx).avgBkg(chIdx, roiIdx, zIdx));

                    end % end ROI loop
                end % End Z loop
                % flip the Z stack if it was acquired from bottom to top of the gel
                if StackResults(stackIdx).flipStack
                    StackResults(stackIdx).Z = flip(StackResults(stackIdx).Z);
                end
                
                % compute the noise as the standard deviation of the
                % background intensity with Z)
                StackResults(stackIdx).bkgNoiseMagnitude(chIdx, roiIdx) = std(StackResults(stackIdx).avgBkg(chIdx, roiIdx, :));

                % use the TopZ value as the relative Z location
                StackResults(stackIdx).Z = StackResults(stackIdx).Z-(double(StackResults(stackIdx).topZ)*double(StackResults(stackIdx).voxelSizeZdefaultValue));
        end % End colour loop

        for roiIdx = 1:StackResults(stackIdx).numROIs
            for chIdx = 1:numChannels
                % Calculate the normalized and minimum-subtracted profile and plot 
                % StackResults(stackIdx).IsumMinSub(chIdx, roiIdx, :) = ...
                %    (StackResults(stackIdx).Isum(chIdx, roiIdx, :)-min(StackResults(stackIdx).Isum(chIdx, roiIdx, :)));
                h = figure(stackIdx); hold on;
                plot(StackResults(stackIdx).Z, squeeze(StackResults(stackIdx).IsumBkgSub(chIdx, roiIdx, :)), cell2mat(defaultColours(chIdx)), 'Linewidth', 2);  
                xlim([min(StackResults(stackIdx).Z), max(StackResults(stackIdx).Z)])
                
                % Save the X-Z summed images
%                 if makeXZPlots
%                     % get the summed intensity values
%                     testImg = sgautoscale(squeeze(StackResults(stackIdx).xzSumImg(chIdx, roiIdx, :, :)));
%                     testSz = size(testImg);
%                     xPlot = [testSz(1).*double(StackResults(stackIdx).voxelSizeXYdefaultValue),0];
%                     zPlot = [testSz(2).*double(StackResults(stackIdx).voxelSizeZdefaultValue),0];
%                     h = figure; imagesc(zPlot, xPlot, testImg, [0, max(testImg(:))/3]); pause(0.5);
%                     title(['stack ', num2str(stackIdx), ', roi ', num2str(roiIdx), ', ', cell2mat(proteins(chIdx))]);
%                     pause(0.5);
%                     pos = get(h, 'Position');
%                     pos(4) = pos(3).*(max(xPlot)/max(zPlot));
%                     set(h, 'Position', pos)
%                     axis image
%                     saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZ.svg']);
%                 end
            end
        end

        % Add legend, x, and y labels
        xlabel('Z depth (\mum)')
        ylabel('Total fluor. intensity (AFU)')
        set(h, 'Position', [150 600 250 140]);

        % Set title
        title(StackResults(stackIdx).title);

        % Save
        saveas(h, [StackResults(stackIdx).outputFilenameBase, 'IProfilesNoLeg.fig']);
        saveas(h, [StackResults(stackIdx).outputFilenameBase, 'IProfilesNoLeg.svg']);
        legend(proteins, 'Location', 'Northeast')
        saveas(h, [StackResults(stackIdx).outputFilenameBase, 'IProfiles.fig']);
        saveas(h, [StackResults(stackIdx).outputFilenameBase, 'IProfiles.svg']);

        % close reader
        close(reader);
    end % End check loop
    
    % make the XZ plots
    for roiIdx = 1:StackResults(stackIdx).numROIs
        %for chIdx = 1:numChannels
        for chIdx = 1:2   
            % Save the X-Z summed images
            if makeXZPlots
                % get the summed intensity values
                testImg = (squeeze(StackResults(stackIdx).xzSumImg(chIdx, roiIdx, :, 1:end-StackResults(stackIdx).topZ)))./5000;
                testImgSub = (squeeze(StackResults(stackIdx).xzSumImgBkgSub(chIdx, roiIdx, :, 1:end-StackResults(stackIdx).topZ)))./5000;

                testSz = size(testImg);
                xPlot = [(testSz(1)-1).*double(StackResults(stackIdx).voxelSizeXYdefaultValue),0];
                zPlot = [(testSz(2)-1).*double(StackResults(stackIdx).voxelSizeZdefaultValue),0];
                x = (xPlot(1):-double(StackResults(stackIdx).voxelSizeXYdefaultValue):xPlot(2))';
                z = zPlot(1):-double(StackResults(stackIdx).voxelSizeZdefaultValue):zPlot(2);
                h = figure; imagesc(zPlot, xPlot, testImg, [0, max(testImg(:))]); pause(0.5);
                title(['stack ', num2str(stackIdx), ', roi ', num2str(roiIdx), ', ', cell2mat(proteins(chIdx))]);
                pause(0.5);
                pos = get(h, 'Position');
                pos(4) = pos(3).*(max(xPlot)/max(zPlot));
                set(h, 'Position', pos)
                axis image
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZ.svg']);
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZ.fig']);
                clf; imagesc(zPlot, xPlot, testImgSub, [0, max(testImgSub(:))]); pause(0.5); axis image;
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZSub.svg']);
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZSub.fig']);
                
                close(h); h = figure; contour3(z, x, testImg, 30); title(['stack ', num2str(stackIdx), ', roi ', num2str(roiIdx), ', ', cell2mat(proteins(chIdx))]); axis image; pause(0.5);
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZContour.svg']);
                clf; contour3(z, x, testImgSub, 30); title(['stack ', num2str(stackIdx), ', roi ', num2str(roiIdx), ', ', cell2mat(proteins(chIdx))]); axis image; pause(0.5);
                saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZContourSub.svg']);

                close(h)
                %saveas(h, [StackResults(stackIdx).outputFilenameBase, '_ROI', num2str(roiIdx), '_', cell2mat(proteins(chIdx)), 'XZContour.fig']);
                %imwrite(testImg, [StackResults(stackIdx).outputFilenameBase, 'ROI', num2str(roiIdx), '_ZX.tif'], 'tif');
                %imwrite(testImg, [StackResults(stackIdx).outputFilenameBase, 'ROI', num2str(roiIdx), '_ZX.png'], 'png');
            end
        end
    end
end % end stack loop
%% 

% Gaussian fit the peaks and extract the migration distances and peak
% widths
for stackIdx = 1:numStacks
    if ~isfield(StackResults(stackIdx), 'intensityFitBounds') || isempty(StackResults(stackIdx).intensityFitBounds) || ~any(any(StackResults(stackIdx).intensityFitBounds > 0))

        % fit the peaks and measure the migration distances and peak widths 
            % for the protein bands
            % set up some of the fields that the measuremigration function needs
            StackResults(stackIdx).stackIdx = stackIdx;
            StackResults(stackIdx).peaknames = proteins;

            % initialize the fields that measuremigration will assign
            StackResults(stackIdx).intPkFitA = 0;
            StackResults(stackIdx).intPkFitMu = 0;
            StackResults(stackIdx).intPkFitSigma = 0;
            StackResults(stackIdx).intPkFitRSq = 0;
            StackResults(stackIdx).intensityFitBounds = 0;
            StackResults(stackIdx).intPkSNR = 0;
            StackResults(stackIdx).intPkSNRReg = 0;

            %StackResults(stackIdx) = measuremigration(StackResults(stackIdx));
            tmp = measuremigration(StackResults(stackIdx));
            StackResults(stackIdx) = tmp;
            save(dataOutputFilename, 'StackResults');
    end
end

% Find the X-Y resolution at the migration distance (Gaussian fit max)

for stackIdx = 1:numStacks
    if ~isfield(StackResults(stackIdx), 'intPkFitXYRes') || isempty(StackResults(stackIdx).intPkFitXYRes) || ~any(StackResults(stackIdx).intPkFitXYRes(:) > 0)
        % fit the peaks and measure the migration distances and peak widths 
            % for the protein bands
            % set up some of the fields that the measuremigration function needs
            StackResults(stackIdx).stackIdx = stackIdx;
            StackResults(stackIdx).peaknames = proteins;

            % initialize the fields that findXYres will assign
            StackResults(stackIdx).intPkFitXRes = 0;
            StackResults(stackIdx).intPkFitYRes = 0;
            StackResults(stackIdx).intPkFitXYRes = 0;
            StackResults(stackIdx).intPkFitXRSq = 0;
            StackResults(stackIdx).intPkFitYRSq = 0;

            %StackResults(stackIdx) = measuremigration(StackResults(stackIdx));
            tmp = findXYres(StackResults(stackIdx));
            StackResults(stackIdx) = tmp;
            save(dataOutputFilename, 'StackResults');
    end
end

% % collect all of the migration distance data for the same experimental
% % conditions
close(reader);
% % find all of the unique experimental conditions from the title field
uniqueConditions = unique({StackResults.title});
numConditions = length(uniqueConditions);
% % loop through the unique experimental conditions, collecting the data for
% % each stack of that condition
for condIdx = 1:numConditions
    % find all the stacks of this experimental condition
    matchingStacks = find(strcmp(uniqueConditions(condIdx), {StackResults.title})>0);
    
    % loop through the matching stacks and collect the data
    for chIdx = 1:numChannels
        for pkIdx = 1:numPeaks
            collectedMu = [];
            collectedSigma = [];
            collectedXY = [];
            collectedSNR = [];
            for stackIdx = 1:5 %stackIdx = matchingStacks
                collectedMu = [collectedMu;nonzeros(squeeze(StackResults(stackIdx).intPkFitMu(chIdx,pkIdx,:)))];
                collectedSigma = [collectedSigma;nonzeros(squeeze(StackResults(stackIdx).intPkFitSigma(chIdx,pkIdx,:)))];
                collectedXY = [collectedXY;nonzeros(squeeze(StackResults(stackIdx).intPkFitXYRes(chIdx,pkIdx,:)))];
                collectedSNR = [collectedSNR;nonzeros(squeeze(StackResults(stackIdx).intPkSNR(chIdx,pkIdx,(StackResults(stackIdx).intPkFitMu(chIdx,pkIdx,:) ~= 0))))]; 
            end
            migrationCollected(condIdx, chIdx) = {collectedMu};
            sigmaCollected(condIdx, chIdx) = {collectedSigma};
            xyResCollected(condIdx, chIdx) = {collectedXY};
            snrCollected(condIdx, chIdx) = {collectedSNR};
            
            migrationMean(condIdx, chIdx) = mean(collectedMu);
            migrationErr(condIdx, chIdx) = std(collectedMu);
            
            sigmaMean(condIdx, chIdx) = mean(collectedSigma);
            collSigmaSq = collectedSigma.^2;
            sigmaSqMean(condIdx, chIdx) = mean(collSigmaSq);
            sigmaErr(condIdx, chIdx) = std(collectedSigma);
            sigmaSqErr(condIdx, chIdx) = std(collSigmaSq);
            
            xyResMean(condIdx, chIdx) = mean(collectedXY);
            xyResSq = collectedXY.^2;
            xyResSqMean(condIdx, chIdx) = mean(xyResSq);
            xyResErr(condIdx, chIdx) = std(collectedXY);
            xyResSqErr(condIdx, chIdx) = std(xyResSq);
            
        end
    end

end
save(dataOutputFilename, 'StackResults', 'migrationCollected', 'sigmaCollected', 'xyResCollected');
figure; plotSpread([migrationCollected(1), migrationCollected(2)], ...
    'categoryLabels', [{'PTB1'},{'GAPDH'}], 'xNames', [{'PTBP1'},{'GAPDH'}], ...
    'categoryColors', [{'b'},{'g'}]);
ylim([0,425]);
ylabel('Migration distance (\mum)');

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
%     figure(hMu); plot(muX, muPlot, ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('EP time (s)'); ylabel('Migration distance (\mum)');
%     xlim([0, 25]);
%     figure(hSigma); plot(sigmaX, sigmaPlot, ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('Diffusion time (s)'); ylabel('\sigma (\mum)');
%     %xlabel('Diffusion time (s)'); ylabel('\sigma^2 (\mum^2)');
%     xlim([0, 35]);
%     figure(hXY); plot(xyX, xyPlot, ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12); 
%     xlabel('Diffusion time (s)'); ylabel('XY peak width (\mum)');
%     %xlabel('Diffusion time (s)'); ylabel('Squared XY peak width (\mum^2)');
%     xlim([0, 35]);
% 
% %     figure(hMu); errorbar(tEP, migrationMean(:, chIdx), migrationErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('EP time (s)'); ylabel('Migration distance (\mum)');
% %     xlim([0, 25]);
% %     figure(hSigma); errorbar(tDiff, sigmaSqMean(:, chIdx), sigmaSqErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('Diffusion time (s)'); ylabel('\sigma^2 (\mum^2)');
% %     xlim([0, 35]);
% %     figure(hXY); errorbar(tDiff, xyResSqMean(:, chIdx), xyResSqErr(:, chIdx), ['.',cell2mat(defaultColours(chIdx))], 'Linestyle', 'none', 'MarkerSize', 12, 'FontSize', 14); 
% %     xlabel('EP time (s)'); ylabel('Squared XY peak width (\mum^2)');
% %     xlim([0, 35]);
% end
% figure(hMu); hold on;
% x = 0:25;
% plot(x, 30*x+1.5, '-.c');
% plot(x, 49*x-210, '-.m');
% plot(x, 24*x-75, '-.r');
% 
