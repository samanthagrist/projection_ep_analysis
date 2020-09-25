%*************************************************************************%
%**           Calculate ROI intensities as a function of Z              **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2018-03-01

function StackResults = calcintensities(StackResults)
% Calculate the summed ROI intensities for all proteins and all tracked 
% ROIs over Z.  Also find the intensity profile for each ROI and fit to a
% Gaussian, extracting the peak width, location, amplitude, and R squared.

proteins = StackResults.proteins;
% set up the bioFormats reader and get the number of planes in the stack
% and number of Z slices
reader = StackResults.imgReader;
reader.setSeries(0);
planes = reader.getImageCount;
totalZ = planes/length(proteins);
numZ = StackResults.numZ;


% calculate the intensity changes
% initialize the output data and video files
resultsEntries = 0;
StackResults.intensityResults = zeros(4, (numZ-StackResults.topZ)*length(proteins)*StackResults.numGoodParticles);

% loop through the colour channels, (assume that the proteins all
    % travelled the same trajectory), perform the ROI calculations for
    % each
for i = 1:length(proteins)
    StackResults.currProtein = i;
    
    % open the video files
    for particleIdx = 1:StackResults.numGoodParticles
        %disp(['beforeVideoWriter1, particleIdx=', num2str(particleIdx)])
        outFile = [StackResults.outputFilenameBase, 'C', num2str(i), '_ROI', num2str(particleIdx), '.avi'];
        %disp(['filename=', outFile])
        vOut(particleIdx) = VideoWriter(outFile);
        open(vOut(particleIdx))
    end
    %disp('beforeVideoWriter2')
    vOutFull = VideoWriter([StackResults.outputFilenameBase, 'C', num2str(i), '.avi']);
    open(vOutFull)
    
       
    % loop through the Z stack
    for j = StackResults.topZ:numZ
        StackResults.currZ = j;
        
        % read in the image
        if StackResults.flipStack
            sliceNum = totalZ-StackResults.currZ + 1;
        else
            sliceNum = StackResults.currZ;
        end
        iPlane = reader.getIndex(sliceNum-1, i-1, 0) + 1;
        StackResults.A =  bfGetPlane(reader, iPlane);
        StackResults.AScaled = uint8(double(double(StackResults.A)*255/double(max(max(StackResults.A)))));

        % loop through the found ROIs
        disp(['looping through ', num2str(StackResults.numGoodParticles), ' particles'])
        for particleIdx = 1:StackResults.numGoodParticles
            resultsEntries = resultsEntries + 1;
            % for this ROI, calculate the intensity change over time and save
            trackSet = cell2mat(StackResults.trackSetArray(particleIdx));
            % initialize the centroid position to the well location
            if j == StackResults.topZ
                particleCentroid = [trackSet(find(trackSet(:,3) == StackResults.topZ-1), 1:2)];
            end
            
            % if we were able to track this centroid all the way to this Z, use
            % the tracked position.  
            if ~isempty(find(trackSet(:,3) == j, 1))
                particleCentroid = [trackSet(find(trackSet(:,3) == j), 1:2)];
            else % Otherwise interpolate between the known tracked positions
                previousPositionIdx = find(trackSet(:,3) < j, 1, 'last');
                previousPosition = trackSet(previousPositionIdx,3);
                nextPositionIdx = find(trackSet(:,3) > j, 1, 'first');
                nextPosition = trackSet(nextPositionIdx,3);
                interpProportion = (j-previousPosition)/(nextPosition-previousPosition);
                previousX = trackSet(previousPositionIdx,1);
                previousY = trackSet(previousPositionIdx,2);
                nextX = trackSet(nextPositionIdx,1);
                nextY = trackSet(nextPositionIdx,2);
                tempX = previousX + interpProportion*(nextX-previousX);
                tempY = previousY + interpProportion*(nextY-previousY); 
                particleCentroid = [tempX, tempY];
            end
            
            % Take the ROI for this Z, this protein, and this protein spot
            currentROI = StackResults.A(uint16(particleCentroid(2)-StackResults.halfROISz):uint16(particleCentroid(2)+StackResults.halfROISz), ...
                uint16(particleCentroid(1)-StackResults.halfROISz):uint16(particleCentroid(1)+StackResults.halfROISz));
            currentROIScaled = StackResults.AScaled(uint16(particleCentroid(2)-StackResults.halfROISz):uint16(particleCentroid(2)+StackResults.halfROISz), ...
                uint16(particleCentroid(1)-StackResults.halfROISz):uint16(particleCentroid(1)+StackResults.halfROISz));
            
            % Write the ROI image to the video file
            writeVideo(vOut(particleIdx),currentROIScaled)
            
            % Sum up the intensity in the ROI and store
            StackResults.intensityResults(1, resultsEntries) = i; % first column is protein
            StackResults.intensityResults(2, resultsEntries) = j; % second column is Z
            StackResults.intensityResults(3, resultsEntries) = particleIdx; % third column is ROI number
            StackResults.intensityResults(4, resultsEntries) = sum(sum(double(currentROI)-StackResults.bkgMean(i))); % last column is summed intensity in ROI
            
            StackResults.summedIntensity(i, j, particleIdx) = sum(sum(double(currentROI)-StackResults.bkgMean(i)));
            StackResults.zLoc(j) = (j-StackResults.topZ)*StackResults.zSpacing;
            
            % if the summed intensity is over the threshold value, take
            % a line profile and Gaussian fit to estimate the width of the 
            % protein spot.
            StackResults.summedIntensity(i, j, particleIdx)
            if (StackResults.summedIntensity(i, j, particleIdx) > StackResults.IThresh) && StackResults.ESTSZ
                for dimIdx = 1:2 % fit the peak width in x and y
                    %y = sum((double(currentROI)-StackResults.bkgMean(i)),dimIdx);
                    y = sum(double(currentROI),dimIdx);
                    y = y - min(y);
                    if ~iscolumn(y)
                        y = y';
                    end
                    x = (1:length(y))'.*StackResults.xyScale;

                    %figure(300); plot(x,y);
                    
                    % fit each ROI cross-section to a Gaussian to find the peak location
                    fit_type = 'gauss1';
                    fit_options = fitoptions(fit_type); 


                    % Set the sigma bounds
                    sigma_min = 0;
                    sigma_max = 2*length(y);

                    % Set the peak center bounds
                    x_min = min(x);
                    x_max = max(x);

                    % Set the ampitude bounds
                    a_min = 0;
                    a_max = max(y);

                    % set the upper and lower bounds. correct for difference in c and
                    % sigma terms
                    fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
                    fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];

                    % Fit the peaks
                    [fit_object, gof] = fit(x, y, fit_type, fit_options);
                    figure(300); plot(fit_object, x, y);
                    title(['Protein ', num2str(i), 'Z ', num2str(j)])

                    % Get the coefficients
                    fit_coeffs = coeffvalues(fit_object);
                    StackResults.xyFitRSq(i,j,particleIdx,dimIdx) = gof.rsquare;
                    gof.rsquare
                    if gof.rsquare > 0.6
                        StackResults.xyFitA(i,j,particleIdx,dimIdx) = fit_coeffs(1);
                        StackResults.xyFitMu(i,j,particleIdx,dimIdx) = fit_coeffs(2);
                        StackResults.xyFitSigma(i,j,particleIdx,dimIdx) = fit_coeffs(3)/sqrt(2);
                        
                    else % we don't have a reliable fit.  Set to zero.
                        StackResults.xyFitA(i,j,particleIdx,dimIdx) = 0;
                        StackResults.xyFitMu(i,j,particleIdx,dimIdx) = 0;
                        StackResults.xyFitSigma(i,j,particleIdx,dimIdx) = 0;
                    end
                    
                end
            else
                StackResults.xyFitA(i,j,particleIdx,1:2) = [0,0];
                StackResults.xyFitMu(i,j,particleIdx,1:2) = [0,0];
                StackResults.xyFitSigma(i,j,particleIdx,1:2) = [0,0];
                StackResults.xyFitRSq(i,j,particleIdx,1:2) = [0,0];
            end
        end
        
        % burn the trajectories into the image and save the output video
        trackCoordsReduced = StackResults.trackCoordsArray(find(~cellfun(@isempty,StackResults.trackCoordsArray)));
        Idisp = cat(3, StackResults.AScaled, StackResults.AScaled, StackResults.AScaled);
        %trackCoordsReduced
        Idisp = insertShape(Idisp, 'Line', trackCoordsReduced, 'Color', 'r');
        %figure(1); imshow(Idisp);
        writeVideo(vOutFull,Idisp)
    end
    
    % close the videos
    for particleIdx = 1:StackResults.numGoodParticles
        close(vOut(particleIdx))
    end
    close(vOutFull)
end
end

