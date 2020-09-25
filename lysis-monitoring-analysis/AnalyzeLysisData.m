%**************************************************************************%
%*                                                                        *%
%*         Analyze fluorescent cell lysis integrated intensities          *%
%* Segments, tracks, and integrates particle intensities in a time series *%
%*                                                                        *%
%*                           Samantha M. Grist                            *%
%*                          v0.1 September 2017                           *%
%*                                                                        *%
%**************************************************************************%

function [tArray, normParticleIntensity, wellIntensity, maxIntensity, normWellIntensity, normMaxIntensity, proteinWidth, meanBkg] = ...
    AnalyzeLysisData(filename, outputFilename, excludeFrames, ...
    tInterval, integrationRadius, wellRadius, dbg, fudge, pixScale, figNum)

% This function requires the tracking library written by Daniel Blair and
% Eric Dufresne, available: http://site.physics.georgetown.edu/matlab/index.html

% INPUTS
% filename: String containing the name of the input AVI data
% outputFilename: String containing the name of the output AVI file
% (modified for .fig file)
% excludeFrames frames to exclude from analysis (focus drift)
% tInterval: time interval, minutes
% integrationRadius: the integration radius, in pixels, for each cell 
% dbg: flag to indicate whether we want to plot debug data
% fudge: thresholding fudge factor.  Decrease to increase sensitivity.
% pixScale: scale of the image (microns/pixel)
% figNum: figure number to which to plot
% OUTPUTS
% tArray: double array containing the time points
% normParticleIntensity: double array containing all of the particle intensities
% for all of the time points

% add the Blair/Dufresne tracking code to the path
addpath('./tracking/')

% initialize the variables
frameNum = 0;
framesAnalyzed = 0;
info = imfinfo(filename);
num_images = numel(info);
vOut = VideoWriter(outputFilename);
open(vOut)

% preallocate for speed
h = info(1).Height;
w = info(1).Width;
tArray = ones(num_images-length(excludeFrames), 1);
segArray = ones(h, w, num_images-length(excludeFrames));
imgArray = ones(h, w, num_images-length(excludeFrames));

skipCtr = 0;
for k = 1:num_images
    frameNum = frameNum+1;
    % read in frame
    I=imread(filename, k, 'Info', info);
    
    % if the frame is not in the exclusion vector, analyze it
    if ~any(excludeFrames == frameNum)
        framesAnalyzed = framesAnalyzed+1;
        disp(['Frame number = ', num2str(frameNum)]);
        
        % update the time variable
        tArray(framesAnalyzed) = frameNum*tInterval;
        
        Ifilt = medfilt2(I); % median filter to remove noise
        % threshold to binarize
        BW = imbinarize(Ifilt,'adaptive','ForegroundPolarity','bright','Sensitivity',fudge);
        %thresh = graythresh(Ifilt);
        %BW = imbinarize(I, thresh*fudge);

        % close up the particles
        se = strel('disk',double(uint8(wellRadius/6)));
        BW = imclose(BW, se);
        
        % remove small noise
        se = strel('disk',double(uint8(wellRadius/4)));
        BW = imopen(BW, se);
        BW = imclose(BW, se);
        figure(figNum); imshow(sgautoscale(I));
        figure(figNum); imshow(BW);
        segArray(:,:,framesAnalyzed) = BW;
        imgArray(:,:,framesAnalyzed) = I;
        
        seNew = strel('disk',double(uint8(2*wellRadius)));
        BWBkg = ~imdilate(BW, seNew);
        meanBkg(k) = mean(I(BWBkg));
        
        % find the centroids
        s = regionprops(BW, 'centroid');
        centroids = cat(1, s.Centroid);
        hold on
        if size(centroids,2) > 1 % if we found centroids
            skipCtr = 0; % reset the skip counter
            plot(centroids(:,1),centroids(:,2), 'b*')
            hold off

            % store the centroids for tracking
            centroids = cat(2, centroids, frameNum*ones(size(centroids, 1), 1));
            if framesAnalyzed == 1
                centroidsArray = centroids;
            else
                centroidsArray = cat(1, centroidsArray, centroids);
            end
        else
            skipCtr = skipCtr+1;
            disp(['No centroids found for ', num2str(skipCtr), ' loops.'])
            % if we haven't found centroids in 10 frames, break the loop
            if skipCtr >= 10
                break;
            end
        end
    else
        disp('skipping')
    end
end
%% 

% track the centroids over time
param = struct('mem', 3, 'good', 10, 'dim', 2 ,'quiet', 0);
tracks = trackModified(centroidsArray, 50, param);
[n, xout] = hist(tracks(:,2),20);
stairs(xout, n, 'color', [0.5 0.5 0.5]);
reducedTracks = tracks;

% remove particles that aren't there from the beginning, and store the
% trajectories
savedFrames = setdiff([1:num_images], excludeFrames);
startVal = min(savedFrames);
numGoodParticles = 0;
for particleIdx = 1:max(tracks(:,4))
    % focus on one particle trajectory  reshape the coordinates for
    % plotting
    trackSet = tracks((tracks(:,4) == particleIdx),:);
    xyPositions = trackSet(:, 1:2);
    minPos = min(xyPositions(:));
    maxPos = max(xyPositions(:));
    
    % if we found this particle in the beginning and it's not too close 
    % to the image edge, store its trajectory
    if (min(trackSet(:,3)) == startVal) && ...
            (minPos > ceil(integrationRadius)) && (maxPos < min([h,w])-ceil(integrationRadius))
        numGoodParticles = numGoodParticles+1;
        trackSetArray(numGoodParticles) = {trackSet};
        trackCoords = [trackSet(:, 1)'; trackSet(:, 2)'];
        trackCoords = trackCoords(:)';
        trackCoordsArray(numGoodParticles) = {trackCoords};
        %trackCoords((1:length(trackSet(:,1)))*2 - 1) = trackSet(:,1);
        %trackCoords((1:length(trackSet(:,2)))*2) = trackSet(:,2);

        % plot command to plot over a figure window instead of burning into
        % the image
        %hold on;
        %plot(tracks(find(tracks(:,4) == particleIdx),1), tracks(find(tracks(:,4) == particleIdx),2));
    else % this particle appeared partway through - ignore it.
        disp(['Removing ', num2str(particleIdx)])
        reducedTracks(reducedTracks(:,4) == particleIdx, :) = [];
    end
end

disp(['Found and tracked ', num2str(numGoodParticles), ' good particles.'])

% calculate the intensity change
particleIntensity = zeros(length(savedFrames), numGoodParticles);
normParticleIntensity = zeros(length(savedFrames), numGoodParticles);
maxIntensity = zeros(length(savedFrames), numGoodParticles);
wellIntensity = zeros(length(savedFrames), numGoodParticles);
normMaxIntensity = zeros(length(savedFrames), numGoodParticles);
normWellIntensity = zeros(length(savedFrames), numGoodParticles);
for particleIdx = 1:numGoodParticles
    disp(['******Particle ', num2str(particleIdx)])
    % for this particle, calculate the intensity change over time and save
    trackSet = reducedTracks((reducedTracks(:,4) == particleIdx),:);
    for frameIdx = min(trackSet(:, 3)'):max(trackSet(:, 3)')
        [tmp, trackIdx] = min(abs(trackSet(:,3) - frameIdx));
        disp(['Frame index = ', num2str(frameIdx), ', track index = ', num2str(trackIdx)])
        particleCentroid = [trackSet(trackIdx, 1:2), integrationRadius];
        particleCentroidWell = [trackSet(trackIdx, 1:2), wellRadius];
        savedFrameIdx = find(savedFrames == frameIdx);
        if dbg
            Idisp = imgArray(:,:,savedFrameIdx)./max(max(imgArray(:,:,savedFrameIdx)));
            Idisp = insertShape(Idisp, 'Circle', particleCentroid);
            figure(figNum); imshow(Idisp);
            pause(0.1);
        end
        intArea = im2bw(insertShape(zeros(h, w), 'FilledCircle', particleCentroid, 'Color', 'w'));
        intAreaWell = im2bw(insertShape(zeros(h, w), 'FilledCircle', particleCentroidWell, 'Color', 'w'));
        lineProfileRegion = im2bw(insertShape(zeros(h, w), 'FilledRectangle', [particleCentroid(1:2)-integrationRadius, 2*integrationRadius, 2*integrationRadius], 'Color', 'w'));
        
        % quantify metrics of interest
        I = imgArray(:,:,savedFrameIdx);
        % calculated the integrated and maximum intensities of interest
        particleIntensity(savedFrameIdx, particleIdx) = sum(double(I(intArea ~= 0)))-meanBkg(frameIdx)*length(nonzeros(intArea ~= 0));
        wellIntensity(savedFrameIdx, particleIdx) = sum(double(I(intAreaWell ~= 0)))-meanBkg(frameIdx)*length(nonzeros(intAreaWell ~= 0));
%         size(I(intArea ~= 0))
%         max(I(intArea ~= 0))
%         max(Ifilt(intArea ~= 0))
%         pause;
        maxIntensity(savedFrameIdx, particleIdx) = max(I(intArea ~= 0))-meanBkg(frameIdx); % find the maximum intensity in the region of the filtered image
        disp(['MaxIntensity = ', num2str(maxIntensity(savedFrameIdx, particleIdx))])
        normParticleIntensity(savedFrameIdx, particleIdx) = double(particleIntensity(savedFrameIdx, particleIdx))/double(particleIntensity((savedFrames == startVal), particleIdx));
        normWellIntensity(savedFrameIdx, particleIdx) = double(wellIntensity(savedFrameIdx, particleIdx))/double(wellIntensity((savedFrames == startVal), particleIdx));
        normMaxIntensity(savedFrameIdx, particleIdx) = double(maxIntensity(savedFrameIdx, particleIdx))/double(maxIntensity((savedFrames == startVal), particleIdx));
        disp(['NormMaxIntensity = ', num2str(normMaxIntensity(savedFrameIdx, particleIdx))])
        % Now calculate the peak widths (via fitting to line profiles)
%         profileRegion = I(lineProfileRegion ~= 0);
%         size(profileRegion)
        maskedI = I.*lineProfileRegion;
        xProfile = sum(maskedI, 2);
        xProfile = nonzeros(xProfile)-min(nonzeros(xProfile));
        yProfile = sum(maskedI, 1);
        yProfile = nonzeros(yProfile)-min(nonzeros(yProfile));
        x = ([1:length(xProfile)].*pixScale)';
%         size(x)
%         size(xProfile)ma
%         size(yProfile)
%         figure; imshow(lineProfileRegion); pause;
        [a, mu, sigma, rsq] = sgGaussFit(x, xProfile);
        [a2, mu2, sigma2, rsq2] = sgGaussFit(x, yProfile);
        disp(['x Rsq = ', num2str(rsq), ', y Rsq = ', num2str(rsq2)]);
        proteinWidth(savedFrameIdx, particleIdx) = mean([sigma, sigma2]);
    end
end

h = figure(figNum); plot(tArray, normParticleIntensity);
xlabel('Time (seconds)')
ylabel('Particle intensity (normalized to initial)')
set(h,'DefaultAxesColorOrder',jet(length(trackCoordsArray)));
saveas(h, strrep(outputFilename, '.avi', '_CellIntensities.fig'), 'fig') 
saveas(h, strrep(outputFilename, '.avi', '_CellIntensities.svg'), 'svg') 

% burn the trajectories into the image and save the output video
for frameIdx = 1:length(tArray)
    Idisp = imgArray(:,:,frameIdx)./max(max(imgArray(:,:,frameIdx)));
    Idisp = insertShape(Idisp, 'Line', trackCoordsArray, 'Color', jet(length(trackCoordsArray)));
    writeVideo(vOut,Idisp)
end
close(vOut)

end