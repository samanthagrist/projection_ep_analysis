%*************************************************************************%
%**    Find the X-Y spot size of each protein at its Z-peak location    **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2019-07-08
% for simple Z profile analysis
 
function [StackResults] = findXYres(StackResults)
colIdx = StackResults.defaultColours;
disp('FIND XY RES!')

proteins = StackResults.proteins;
numZ = StackResults.numZ;
numGoodParticles = StackResults.numROIs;
roiSz = StackResults.roiSz;
wBkgXY = StackResults.wBkgXY;

% initialize the reader
reader = bfGetReader(StackResults.inputFilename); % use the bioformats toolbox to create the image reader for this stack
reader.setSeries(0);
planes = reader.getImageCount;

% loop through the colour channels
for colourIdx = 1:length(proteins)
    % loop through all the peaks defined for this colour channel
    for pkIdx = 1:size(StackResults.peaknames, 2)
        
        if ~isempty(cell2mat(StackResults.peaknames(colourIdx,pkIdx))) % if the name of this peak is not empty run the calc
            % fit each peak to a Gaussian for each ROI
            for roiIdx = 1:numGoodParticles
                disp(['Colour ', num2str(colourIdx), ', peak ', num2str(pkIdx), ', ROI ', num2str(roiIdx)]);
                % get the peak location for this peak in this colour channel, from
                % the data generated by measuremigration
                x = StackResults.Z';
                pkLoc = StackResults.intPkFitMu(colourIdx,pkIdx,roiIdx);
                [tmp, pkLocI] = min(abs(double(StackResults.Z) - double(pkLoc)));
                
                % load in the Z-slice image corresponding to the migration peak
                disp(['Ch. ', num2str(colourIdx), ', Z ', num2str(pkLocI)])
                % Read in the image
                iPlane = reader.getIndex(pkLocI-1, colourIdx-1, 0) + 1;
                A = bfGetPlane(reader, iPlane);
                ASz = size(A);
                
                % isolate the ROI
                [xStart, yStart, xEnd, yEnd] = findROILims(StackResults, roiIdx, roiSz, ASz, 1);
                disp([num2str(xStart), ' ', num2str(xEnd), ' ', num2str(yStart), ' ', num2str(yEnd)])
          

                Aregion = A(xStart:xEnd, yStart:yEnd);
                
                % estimate the x- and y- directional background
                %w = gausswin(10)/sum(gausswin(10));
                gaussFiltSz = 10;
                XBkgTop = double(Aregion(:, 1:wBkgXY));
                XBkgBottom = double(Aregion(:, (roiSz+1-wBkgXY):roiSz));
                XBkg = [XBkgTop,XBkgTop];
                XBkgProfNoisy = mean(double(XBkg), 2);
                XBkgProf = mean(double(imgaussfilt(XBkg, gaussFiltSz)), 2);
                YBkgTop = double(Aregion(1:wBkgXY, :));
                YBkgBottom = double(Aregion((roiSz+1-wBkgXY):roiSz, :));
                YBkg = [YBkgTop;YBkgTop];
                YBkgProfNoisy = mean(double(YBkg), 1);
                YBkgProf = mean(double(imgaussfilt(YBkg, gaussFiltSz)), 1);
                
                %figure(500); clf; plot(YBkgProfNoisy, 'b'); hold on; plot(YBkgProf, 'r');
                           
                % create the bkg-subtracted x- and y-averaged profiles
                XAvgProf = mean(double(Aregion), 2)-XBkgProf;
                YAvgProf = mean(double(Aregion), 1)-YBkgProf;
                pixelSz = double(StackResults.voxelSizeXYdefaultValue);
                xyPosIdx = 0:1:(length(XAvgProf)-1);
                relXYPos = xyPosIdx.*pixelSz;
                
                % set up and run the fits
                % fit each ROI cross-section to a Gaussian to find the X-Y peak width
                
                % first find X peak width
                y1 = XAvgProf;
                x = relXYPos';
                %length(x)
                [fit_object_x, gof_x] = gaussFitXY(x,y1);
                
                % next find Y peak width
                y2 = YAvgProf';
                [fit_object_y, gof_y] = gaussFitXY(x,y2);
                
                % save the peak widths, and goodness of fit
                fit_coeffsX = coeffvalues(fit_object_x);
                fit_coeffsY = coeffvalues(fit_object_y);
                sigmaX = fit_coeffsX(3)/sqrt(2);
                sigmaY = fit_coeffsY(3)/sqrt(2);
                
                rsqX = gof_x.rsquare;
                rsqY = gof_y.rsquare;
                
                % plot the fits
                figure(300); subplot(2,2,1, 'replace');
                plot(fit_object_x, x, y1); 
                xlabel('X distance'); ylabel('Average intensity'); title(['X Rsq = ', num2str(rsqX)]);
                subplot(2,2,2, 'replace'); plot(fit_object_y, x, y2); 
                xlabel('Y distance'); ylabel('Average intensity'); title(['Y Rsq = ', num2str(rsqY)]);
                subplot(2,2,3, 'replace'); imshow(sgautoscale(Aregion));
                subplot(2,2,4, 'replace');
                pause(0.01); plot(YBkgProfNoisy, 'b'); hold on; plot(YBkgProf, 'r');

                
                
                % if we only have rsq>0.6 (reliable fit) for one of the directions, use
                % that direction as the peak width.  If we don't have it for either, set peak width to 0.
                if rsqX > 0.6
                    StackResults.intPkFitXRes(colourIdx,pkIdx,roiIdx) = sigmaX;
                    StackResults.intPkFitXRSq(colourIdx,pkIdx,roiIdx) = rsqX;
                    if rsqY > 0.6                  
                        StackResults.intPkFitYRes(colourIdx,pkIdx,roiIdx) = sigmaY;
                        StackResults.intPkFitXYRes(colourIdx,pkIdx,roiIdx) = mean([sigmaX,sigmaY]);
                        StackResults.intPkFitYRSq(colourIdx,pkIdx,roiIdx) = rsqY;
                    else
                        StackResults.intPkFitYRes(colourIdx,pkIdx,roiIdx) = 0;
                        StackResults.intPkFitXYRes(colourIdx,pkIdx,roiIdx) = sigmaX;
                        StackResults.intPkFitYRSq(colourIdx,pkIdx,roiIdx) = rsqY;
                    end
                else
                    StackResults.intPkFitXRes(colourIdx,pkIdx,roiIdx) = 0;
                    StackResults.intPkFitXRSq(colourIdx,pkIdx,roiIdx) = rsqX;
                    if rsqY > 0.6
                        StackResults.intPkFitYRes(colourIdx,pkIdx,roiIdx) = sigmaY;
                        StackResults.intPkFitXYRes(colourIdx,pkIdx,roiIdx) = sigmaY;
                        StackResults.intPkFitYRSq(colourIdx,pkIdx,roiIdx) = rsqY;
                    else
                        StackResults.intPkFitYRes(colourIdx,pkIdx,roiIdx) = 0;
                        StackResults.intPkFitXYRes(colourIdx,pkIdx,roiIdx) = 0;
                        StackResults.intPkFitYRSq(colourIdx,pkIdx,roiIdx) = rsqY;
                    end
                end
                        

%                 % save all information in the master data structure
%                 % Format: (1) Gel number	(2) Gel %T	(3) Gel %R	(4) EP time	(5) EP field ..
%                 % (6) Stack index	(7) ROI index   (8) Protein channel	(9) Peak name ..
%                 % (10) Peak A	(11) Peak Mu	(12) Peak Sigma	(13) Peak Rsq	(14) x res at peak ..
%                 % (15) y res at peak	(16) diffusion time 
%                 CombinedResults.IntensityPeakEntries = CombinedResults.IntensityPeakEntries+1;
%                 if isfield(StackResults, 'tDiff')
%                     dataEntry = {StackResults.gelNumber, StackResults.gelT, StackResults.gelR, ...
%                         StackResults.epTime, StackResults.epField, StackResults.stackIdx, particleIdx, ...
%                         proteins(colourIdx), StackResults.peaknames(colourIdx, pkIdx), a, ...
%                         muDepth, sigmaDepth, rsq, xyres(1), xyres(2), StackResults.tDiff}';
%                 else % backwards compatibility with older sets of data
%                     dataEntry = {StackResults.gelNumber, StackResults.gelT, StackResults.gelR, ...
%                         StackResults.epTime, StackResults.epField, StackResults.stackIdx, particleIdx, ...
%                         proteins(colourIdx), StackResults.peaknames(colourIdx, pkIdx), a, ...
%                         muDepth, sigmaDepth, rsq, xyres(1), xyres(2), StackResults.epTime}';
%                 end
%                 CombinedResults.IntensityPeaks(CombinedResults.IntensityPeakEntries, :) = dataEntry;
            end
                
        end
    end

end
close(reader)
end