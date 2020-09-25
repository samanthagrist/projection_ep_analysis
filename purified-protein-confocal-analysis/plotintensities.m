%*************************************************************************%
%**              Plot ROI intensities as a function of Z                **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2018-03-02

function StackResults = plotintensities(StackResults)
% Plot the summed ROI intensities for all protein channels and all tracked 
% ROIs over Z.  Plot both the averaged signal from all ROIs and the overlay
% of the ROIs, in both raw (background subtracted) intensity and SNR.

colIdx = StackResults.colIdx;

proteins = StackResults.proteins;
numZ = StackResults.numZ;
numGoodParticles = StackResults.numGoodParticles;

% Plot the intensities
for i=1:length(proteins)
    temp = mean(StackResults.summedIntensity(i,:,:), 3);
    StackResults.avgIntensity(i, 1:length(temp)) = temp;
    StackResults.stdIntensity(i, 1:length(temp)) = (std(squeeze(StackResults.summedIntensity(i,:,:))'))';
    
    % Intensity average 
    figure(10); hold on;
    errorbar(StackResults.zLoc(StackResults.topZ:numZ), ...
        StackResults.avgIntensity(i, StackResults.topZ:numZ), ...
        StackResults.stdIntensity(i, StackResults.topZ:numZ), colIdx(i));
    
    % Intensity overlay
    figure(11); hold on;
    for particleIdx = 1:numGoodParticles
        plot(StackResults.zLoc(StackResults.topZ:numZ), StackResults.summedIntensity(i,StackResults.topZ:numZ,particleIdx), colIdx(i))
    end
    
    % SNR average
    % get the size of the region that we summed over
    regLength = 2*StackResults.halfROISz+1;
    regSz = numel(zeros(regLength, regLength));
    figure(12); hold on;
    errorbar(StackResults.zLoc(StackResults.topZ:numZ), ...
        StackResults.avgIntensity(i, StackResults.topZ:numZ)./((regSz)*StackResults.bkgStd(i)), ...
        StackResults.stdIntensity(i, StackResults.topZ:numZ)./((regSz)*StackResults.bkgStd(i)), colIdx(i));
    
    % SNR overlay
    figure(13); hold on;
    for particleIdx = 1:numGoodParticles
        plot(StackResults.zLoc(StackResults.topZ:numZ), ...
            StackResults.summedIntensity(i,StackResults.topZ:numZ,particleIdx)./((regSz)*StackResults.bkgStd(i)), colIdx(i))
    end
    
    % Protein spot peak width
    if StackResults.ESTSZ
        % xyFitMu
        size(StackResults.xyFitMu(i,:,:,:))
        % find the elements where both the x- and y- peak widths are
        % nonzero
        [row1, col1] = find(squeeze(StackResults.xyFitSigma(i,:,:,1)) ~= 0);
        [row2, col2] = find(squeeze(StackResults.xyFitSigma(i,:,:,2)) ~= 0);
        
%          common = find((StackResults.xyFitSigma(i,:,:,2)) ~= 0 & (StackResults.xyFitSigma(i,:,:,1)) ~= 0);
%          size(common)
%          common
%         pause;
%         
        if isrow(row1)
            A1 = [row1', col1'];
            A2 = [row2', col2'];
        else
            A1 = [row1, col1];
            A2 = [row2, col2];
        end
        
%         size(A1)
%         size(A2)
%         nnz(StackResults.xyFitSigma(i,:,:,1))
%         nnz(StackResults.xyFitSigma(i,:,:,2))
        %A1
        %A2
        %i
        if ~isempty(A1) && ~isempty(A2)
        %if ~isempty(common)
            common = intersect(A1, A2, 'rows');

            % only look at the nonzero elements for calculation and plot
            common
            if StackResults.numGoodParticles == 1
                xyCondensed = StackResults.xyFitSigma(i,common(:,2),1,:);
            else
                xyCondensed = StackResults.xyFitSigma(i,common(:,1),common(:,2),:);
            end

            temp = mean(xyCondensed(1,:,:,:), 4);

            StackResults.avgWidth(i, 1:size(temp,2)) = mean(temp, 3);

            StackResults.stdWidth(i, 1:size(temp,2)) = (std(temp, 0, 3));

            figure(14); hold on;
            errorbar(StackResults.zLoc(common(:,1)), ...
                StackResults.avgWidth(i, 1:size(temp,2)), ...
                StackResults.stdWidth(i, 1:size(temp,2)), colIdx(i));
        else
            figure(14); hold on;
            plot(StackResults.zLoc, ...
                zeros(length(StackResults.zLoc), 1), ...
                colIdx(i));
        end
end
h = figure(10);
xlabel('Z depth (\mum)')
ylabel('Summed ROI intensity')
legend(proteins);
set(h,'DefaultAxesColorOrder',jet(length(StackResults.trackCoordsArray)));
saveas(h, [StackResults.outputFilenameBase, '_ROIIntensities.fig'], 'fig') 
saveas(h, [StackResults.outputFilenameBase, '_ROIIntensities.svg'], 'svg')
hold off;

h = figure(11);
xlabel('Z depth (\mum)')
ylabel('Summed ROI intensity')
legend('off');
set(h, 'Position', [500 500 180 180]);
set(h,'DefaultAxesColorOrder',jet(length(StackResults.trackCoordsArray)));
saveas(h, [StackResults.outputFilenameBase, '_ROIIntensitiesAll.fig'], 'fig') 
saveas(h, [StackResults.outputFilenameBase, '_ROIIntensitiesAll.svg'], 'svg')
hold off;

h = figure(12);
xlabel('Z depth (\mum)')
ylabel('Summed ROI SNR')
legend(proteins);
set(h,'DefaultAxesColorOrder',jet(length(StackResults.trackCoordsArray)));
saveas(h, [StackResults.outputFilenameBase, '_ROISNR.fig'], 'fig') 
saveas(h, [StackResults.outputFilenameBase, '_ROISNR.svg'], 'svg')
hold off;

h = figure(13);
xlabel('Z depth (\mum)')
ylabel('Summed ROI SNR')
legend('off');
set(h, 'Position', [500 500 180 180]);
set(h,'DefaultAxesColorOrder',jet(length(StackResults.trackCoordsArray)));
saveas(h, [StackResults.outputFilenameBase, '_ROISNRAll.fig'], 'fig') 
saveas(h, [StackResults.outputFilenameBase, '_ROISNRAll.svg'], 'svg')
hold off;

if StackResults.ESTSZ
    h = figure(14);
    xlabel('Z depth (\mum)')
    ylabel('Average protein spot peak width (\mum)')
    legend(proteins);
    set(h,'DefaultAxesColorOrder',jet(length(StackResults.trackCoordsArray)));
    saveas(h, [StackResults.outputFilenameBase, '_ROIPkWidth.fig'], 'fig') 
    saveas(h, [StackResults.outputFilenameBase, '_ROIPkWidth.svg'], 'svg')
    hold off;
end

end