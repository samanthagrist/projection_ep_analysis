%*************************************************************************%
%**              Track protein spot centroids through Z                 **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2018-03-01

function StackResults = trackcentroids(StackResults)
% find the centroids of the ROI regions in the summed image from the three 
% colour channels, and track their location through Z.  Filter out
% centroids that disappear or move out of the field of view before 3/4 
% through the stack

proteins = StackResults.proteins;
% set up the bioFormats reader and get the number of planes in the stack
% and number of Z slices
reader = StackResults.imgReader;
reader.setSeries(0);
planes = reader.getImageCount;
totalZ = planes/length(proteins);
numZ = StackResults.numZ;


% read in an example image
% set up the plane index from the Z val and colour
iPlane = reader.getIndex(0, 0, 0) + 1;
% read in an example image
I = bfGetPlane(reader, iPlane);
minZ = StackResults.topZ;
maxZ = numZ;
%minZ
%maxZ

% initialize the array of centroids to the starting well locations
if isfield(StackResults, 'wellLocs') && ~isempty(StackResults.wellLocs)
    for proteinIdx = 1:length(proteins)
        wells = StackResults.wellLocs;
        StackResults.centroids(StackResults.exampleProtein,  StackResults.topZ-1, :) = wells;
        centroidsArray = cat(2, cell2mat(wells), (minZ-1)*ones(size(cell2mat(wells), 1), 1));
    end
else
    centroidsArray = [];
end

% for summed image of 3 channels go through Z and segment to find ROIs.  
for j = minZ:maxZ
    StackResults.currZ = j;
    A = uint16(zeros(size(I)));
    % sum up the colour channels to create overlay image
    for proteinIdx = 1:length(proteins)
        StackResults.currProtein = StackResults.exampleProtein;
        % calculate the slice number depending on whether we are flipping
        % the stack
        if StackResults.flipStack
            sliceNum = totalZ-j + 1;
        else
            sliceNum = j;
        end
        %sliceNum
        iPlane = reader.getIndex(sliceNum-1, proteinIdx-1, 0) + 1;
        A = A + bfGetPlane(reader, iPlane);
    end
    
    StackResults.A = A;
    StackResults.AScaled = uint8(double(double(StackResults.A)*255/double(max(max(StackResults.A)))));
    StackResults = segmentImage(StackResults);
    centroids = cell2mat(StackResults.centroids(StackResults.exampleProtein, StackResults.currZ));
  
    % store the centroids for tracking
    centroidsForTracking = cat(2, centroids, j*ones(size(centroids, 1), 1));
%     if j == StackResults.topZ
%         startCentroids
%         centroidsArray = centroidsForTracking;
%     else
    if ~isempty(centroids)
        centroidsArray = cat(1, centroidsArray, centroidsForTracking);
    end
%     end
end


% Run tracking algorithm to find ROI displacements
% track the centroids over time
param = struct('mem', 70, 'good', 10, 'dim', 2 ,'quiet', 0);
tracks = trackModified(centroidsArray, 50, param);
[n, xout] = hist(tracks(:,2),20);
stairs(xout, n, 'color', [0.5 0.5 0.5]);
reducedTracks = tracks;

% remove ROIs that aren't there from the beginning, and store the
% trajectories
StackResults.imgSz = size(StackResults.A);
startVal = StackResults.topZ;
endVal = numZ;
StackResults.numGoodParticles = 0;
numTracks = max(tracks(:,4));
disp(['Found ', num2str(numTracks), ' particles to filter.'])
for particleIdx = 1:numTracks
    % focus on one particle trajectory  reshape the coordinates for
    % plotting
    trackSet = tracks((tracks(:,4) == particleIdx),:);
    
    % Criteria for keeping a track and storing its trajectory: 
        % 1. we found this ROI in the beginning (it was a well)
        % 2. it runs >3/4 of the way through the stack
        % 3. its entire ROI size remains in the field of view throughout the experiment 
        particleIdx
        minIdx = min(trackSet(:,3));
        %minIdx
%         if minIdx == startVal-1 
%             maxIdx = max(trackSet(:,3));
%             maxIdx
%             
%             if maxIdx > 0.8*endVal
%                 minLoc = min(min((trackSet(:, 1:2))))
%                 maxLoc = max(trackSet(:, 1))
%                 size(trackSet)
%             end
%         end
        
    if ((min(trackSet(:,3)) == startVal-1) || (min(trackSet(:,3)) == startVal)) && (max(trackSet(:,3)) > 0.7*endVal)
        if  (min(min((trackSet(:, 1:2)))) > StackResults.halfROISz) && ...
            (max(trackSet(:, 1)) < (StackResults.imgSz(1)-StackResults.halfROISz)) && ...
            (max(trackSet(:, 2)) < (StackResults.imgSz(2)-StackResults.halfROISz))
            StackResults.numGoodParticles = StackResults.numGoodParticles+1;
            StackResults.trackSetArray(StackResults.numGoodParticles) = {trackSet};
            trackCoords = [trackSet(:, 1)'; trackSet(:, 2)'];
            trackCoords = trackCoords(:)';
            StackResults.trackCoordsArray(particleIdx) = {trackCoords};
        else
            disp(['Removing ', num2str(particleIdx), ' (too close to edge)'])
            reducedTracks(reducedTracks(:,4) == particleIdx, :) = [];
        end
    else % this particle appeared partway through - ignore it.
        disp(['Removing ', num2str(particleIdx), ' with topVal ', num2str(min(trackSet(:,3))), ' and bottomVal ', num2str(max(trackSet(:,3)))])
        reducedTracks(reducedTracks(:,4) == particleIdx, :) = [];
    end
end
disp(['Found ', num2str(StackResults.numGoodParticles), ' good particles.'])

end