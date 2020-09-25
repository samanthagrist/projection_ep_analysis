%*************************************************************************%
%**                   Measure migration information                     **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2019-07-08
% for simple Z profile analysis
 
function [StackResults] = measuremigrationLS(StackResults)
% Prompt the user to define the peak boundaries for each peak of interest
% in each colour channel, then fit each peak to a Gaussian for each ROI,
% extracting the  amplitude, mean, peak width, and goodness of fit
% information. At the peak location for each peak, find the x-y resolution 
% peak width.  Save this data in the StackResults object as well as in a 
% master data structure.
colIdx = StackResults.defaultColours;
disp('MEASURE MIGRATION!')
dbg = 0;

proteins = StackResults.proteins;
numZ = StackResults.numZ;
numGoodParticles = StackResults.numROIs;

% loop through the colour channels
for colourIdx = 1:length(proteins)
    % plot the intensity profiles for all ROIs and prompt the user to define
    % the peak boundaries for each peak
    figure(11); clf; hold on;
    for particleIdx = 1:numGoodParticles
        x = squeeze(StackResults.ZScaled(particleIdx, :))';
        y = squeeze(StackResults.IsumBkgSub(colourIdx, particleIdx, :)); 
        %disp([x y])

        plot(x,y,'Color',(colIdx(colourIdx,:)))
    end
    
    % loop through all the peaks defined for this colour channel
    for pkIdx = 1:size(StackResults.peaknames, 2)
        if ~isempty(cell2mat(StackResults(1).peaknames(colourIdx,pkIdx))) % if the name of this peak is not empty run the calc
            figure(11);
            % let the user choose the peak region of interest to fit the Gaussian within
            % Get the limits of the plot
            y_lim = get(gca, 'YLim');
            x_lim = get(gca, 'XLim');
            % prompt the user to pick the peak boundaries
            title(['Stack ', num2str(StackResults.stackIdx), ': click on the left boundary of the ', ...
                cell2mat(StackResults.peaknames(colourIdx, pkIdx)), ' peak.']);
            [x,y,button] = ginput(1);
            leftBound = uint16(x);
            if leftBound < min(x_lim)+1
                leftBound = min(x_lim)+1;
            end
            StackResults.intensityFitBounds(colourIdx, pkIdx, 1) = leftBound;
            % Draw the selected peak boundary
            l1 = line([x, x], y_lim, 'Color', [1, 0, 0]);

            title(['Stack ', num2str(StackResults.stackIdx), ': click on the right boundary of the ', ...
                cell2mat(StackResults.peaknames(colourIdx, pkIdx)), ' peak.']);
            [x,y,button] = ginput(1);
            rightBound = uint16(x);
            if rightBound > max(x_lim)-1
                rightBound = max(x_lim)-1;
            end
            StackResults.intensityFitBounds(colourIdx, pkIdx, 2) = rightBound;
            % Draw the selected peak boundary
            l2 = line([x, x], y_lim, 'Color', [1, 0, 0]);
            leftBound 
            rightBound
            
            % fit each peak to a Gaussian for each ROI
            for particleIdx = 1:numGoodParticles
                % set up and run the fit
                % fit each ROI cross-section to a Gaussian to find the peak location
           
                
                %y = squeeze(StackResults.IsumBkgSub(colourIdx, particleIdx, ...
                    %StackResults.topZ+leftBound:StackResults.topZ+rightBound, particleIdx));
                [tmp,leftBoundIdx] = min(abs(double(StackResults.ZScaled(particleIdx, :)) - double(leftBound)));
                [tmp,rightBoundIdx] = min(abs(double(StackResults.ZScaled(particleIdx, :)) - double(rightBound)));
                %disp([leftBoundIdx, rightBoundIdx])
                y = squeeze(StackResults.IsumBkgSub(colourIdx, particleIdx, rightBoundIdx:leftBoundIdx));
                x = squeeze(StackResults.ZScaled(particleIdx, rightBoundIdx:leftBoundIdx))';
                bkg = squeeze(StackResults.avgBkg(colourIdx, particleIdx, rightBoundIdx:leftBoundIdx));
                bkgNoise = squeeze(StackResults.avgBkgNoiseRegionBkgSub(colourIdx, particleIdx, rightBoundIdx:leftBoundIdx));
                bkgNoiseStd = double(std(bkg))*((StackResults.roiSz+1)^2);
                bkgNoiseStdNoiseRegion = double(std(bkgNoise))*((StackResults.roiSz+1)^2);
                
                
                %length(x)
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
                if max(y) > 0
                    a_max = 1.5*max(y);
                else % if the whole range is below zero, artificially shift it up so we can fit.
                    y = y - min(y);
                    a_max = 1.5*max(y);
                end

                % set the upper and lower bounds. correct for difference in c and
                % sigma terms
                fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
                fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];
%                 fit_options.Lower
%                 fit_options.Upper

                % Fit the peaks
                [fit_object, gof] = fit(x, y, fit_type, fit_options);
                figure(300); clf; plot(fit_object, x, y); pause(0.01);

                % save the amplitude, mean, peak width, and goodness of fit
                fit_coeffs = coeffvalues(fit_object);
                a = fit_coeffs(1);
                mu = fit_coeffs(2);
                %muDepth = double(StackResults.voxelSizeZdefaultValue)*(double(leftBound)+double(mu));
                sigma = fit_coeffs(3)/sqrt(2);
                %sigmaDepth = double(StackResults.voxelSizeZdefaultValue)*double(sigma);
                rsq = gof.rsquare;
                
                % compute the SNR as the peak amplitude divided by
                % the standard deviation of the background Z-profile in the
                % peak region - https://aem.asm.org/content/74/10/2957 Empirical Evaluation of a New Method for Calculating Signal-to-Noise Ratio for Microarray Data Analysis
                % Zhili He, Jizhong Zhou
                % Applied and Environmental Microbiology May 2008, 74 (10) 2957-2966; DOI: 10.1128/AEM.02536-07
                snr_est = a/(2*bkgNoiseStd);
                snr_est2 = a/(2*bkgNoiseStdNoiseRegion);  
                disp(['Rsq = ', num2str(rsq)]);
                disp(['SNR = ', num2str(snr_est2)]);
                if dbg
                    disp(['noise = ', num2str(bkgNoiseStd), ', amplitude = ', num2str(a), ', SNR = ', num2str(snr_est), ', Noise region SNR = ', num2str(snr_est2)]);
                    figure(400); clf; plot(fit_object, x, y); hold on; plot(x, double((bkg-mean(bkg)))*((StackResults.roiSz+1)^2)); 
                    figure(401); clf; plot(fit_object, x, y); hold on; plot(x, double((bkgNoise))*((StackResults.roiSz+1)^2)); pause;
                end
                
               
                if rsq > StackResults.rsqThresh && snr_est2 > 3
                    StackResults.intPkFitA(colourIdx,pkIdx,particleIdx) = a;
                    StackResults.intPkFitMu(colourIdx,pkIdx,particleIdx) = mu;
                    StackResults.intPkFitSigma(colourIdx,pkIdx,particleIdx) = sigma;
                    StackResults.intPkFitRSq(colourIdx,pkIdx,particleIdx) = rsq;
                    StackResults.intPkSNR(colourIdx,pkIdx,particleIdx) = snr_est2;
                    StackResults.intPkSNRReg(colourIdx,pkIdx,particleIdx) = snr_est;
                    

%                     % find the x-y resolution peak width at the peak location
%                     if StackResults.ESTSZ
%                         xyres = ...
%                             StackResults.xyFitSigma(colourIdx, ...
%                             StackResults.topZ+leftBound+uint16(mu),particleIdx,1:2);
% 
%                         StackResults.intPkFitXYRes(colourIdx,pkIdx,particleIdx, 1:2) = xyres;
%                     end
                else
                    StackResults.intPkFitA(colourIdx,pkIdx,particleIdx) = 0;
                    StackResults.intPkFitMu(colourIdx,pkIdx,particleIdx) = 0;
                    StackResults.intPkFitSigma(colourIdx,pkIdx,particleIdx) = 0;
                    StackResults.intPkFitRSq(colourIdx,pkIdx,particleIdx) = rsq;
                    StackResults.intPkSNR(colourIdx,pkIdx,particleIdx) = snr_est2;
                    StackResults.intPkSNRReg(colourIdx,pkIdx,particleIdx) = snr_est;
%                     if StackResults.ESTSZ
%                         StackResults.intPkFitXYRes(colourIdx,pkIdx,particleIdx, 1:2) = 0;
%                     end
                    xyres(1) = 0;
                    xyres(2) = 0;
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
end