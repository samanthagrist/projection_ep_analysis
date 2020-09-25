
%*************************************************************************%
%**                     Plot migration information                      **%
%*************************************************************************%
% Samantha M. Grist
% v.0.1 - 2018-03-06
 
function [StackResults, CombinedResults] = plotmigration(StackResults, CombinedResults)
% Plot the migration information that we collected from the intensity
% stacks
colIdx = StackResults.colIdx;
disp('PLOT MIGRATION!')
figSzX = 520;
figSzY = 260;


proteins = StackResults.proteins;
numZ = StackResults.numZ;
numGoodParticles = StackResults.numGoodParticles;
colIdx = 0;
hMu = [];
hSigma = [];
hXYRes = [];
hMuFit = [];
hSigmaFit = [];
hXYResFit = [];

allPeaks = CombinedResults.IntensityPeaks(:,9);
allPeaks = [allPeaks{:}];
peaks = unique(allPeaks);
for iPeak = 1:length(peaks)
    LegendStrFitMu(iPeak, :) = {[]};
    LegendStrFitSigma(iPeak, :) = {[]};
    LegendStrFitXY(iPeak, :) = {[]};
end


% crazy loop of death to go through the data from all ROIs and images of the same electrophoresis conditions
% ugh I am sure there is a better way of doing this.
gelTs = unique(cell2mat(CombinedResults.IntensityPeaks(:,2)));
% go through all of the gel %Ts
for iT = 1:length(gelTs)
    indT = find(cell2mat(CombinedResults.IntensityPeaks(:,2)) == gelTs(iT)); % subset of the indices with this %T
    gelRs = unique(cell2mat(CombinedResults.IntensityPeaks(indT,3)));
    % go through all the gel %Rs
    for iR = 1:length(gelRs)
        indR = indT(find(cell2mat(CombinedResults.IntensityPeaks(indT,3)) == gelRs(iR))); % subset of the indices with this %rhinohide
        eFields = unique(cell2mat(CombinedResults.IntensityPeaks(indR,5))); 
        % go through all the field strengths
        for iField = 1:length(eFields)
            indE = indR(find(cell2mat(CombinedResults.IntensityPeaks(indR,5)) == eFields(iField))); % subset of the indices with this field
            gelReps = unique(cell2mat(CombinedResults.IntensityPeaks(indE,1)));
            gelReps
            length(gelReps)
            % go through all of the gel replicates
            for igelRep = 1:length(gelReps)
                indRep = indE(find(cell2mat(CombinedResults.IntensityPeaks(indR,1)) == gelReps(igelRep))); % subset of the indices with this gel replicate #
                % set up the plotting colormap and legend string
                if length(gelTs)*length(gelRs)*length(eFields)*length(gelReps)>1
                    colArray = [{[0,0,0]}; {[0.7,0.7,0.7]}; num2cell(jet(length(gelTs)*length(gelRs)*length(eFields)*length(gelReps)+1),2)];
                    suppressLegend = 0;
                else
                    colArray = [{'k'}];
                    suppressLegend = 1;
                end
                colIdx = colIdx+1;
                legendStr(colIdx) = {[num2str(gelTs(iT)), '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm, Replicate ', num2str(gelReps(igelRep))]};
                allPeaks = CombinedResults.IntensityPeaks(indRep,9);
                allPeaks = [allPeaks{:}];
                peaks = unique(allPeaks); 


                % for each protein peak, go through all of the data
                for iPeak = 1:length(peaks)
                    indPk = indRep(find(ismember(allPeaks,peaks(iPeak)))); % subset of the indices with this peak name
                    epTimes = unique(cell2mat(CombinedResults.IntensityPeaks(indPk,4)));

                    % initialize our arrays for the boxplots
                    muArray = [];
                    sigmaArray = [];
                    xyresArray = [];
                    diffTimeArrayXY = []; % use times as a grouping variable.
                    epTimeArrayMu = []; % use times as a grouping variable.
                    diffTimeArraySigma = []; % use times as a grouping variable.
                    
                    % narrow the set of EP times to those that we said we
                    % wanted to plot
                    if isfield(StackResults, 'plotEPtimes')
                        if StackResults.plotEPtimes ~= 0
                            epTimes = intersect(epTimes, StackResults.plotEPtimes);
                        end
                    end

                    % go through all the EP times.  This will be what we plot
                    % against.
                    for iTime = 1:length(epTimes)
                        indTime = indPk(find(cell2mat(CombinedResults.IntensityPeaks(indPk,4)) == epTimes(iTime)));

                        % take all of the ROIs from each of the images and each
                        % of the gels with this condition, and store them to
                        % make the boxplot.
                        tmpMu = {nonzeros(cell2mat(CombinedResults.IntensityPeaks(indTime,11)))};
                        mu(iT, iR, iField, iPeak, iTime, :) = tmpMu;
                        muArray = [muArray; cell2mat(tmpMu)];

                        tmpSigma = {nonzeros(cell2mat(CombinedResults.IntensityPeaks(indTime,12)))};
                        sigNonZInd = find(cell2mat(CombinedResults.IntensityPeaks(indTime,12)) ~= 0);
                        sigma(iT, iR, iField, iPeak, iTime, :) = tmpSigma;
                        sigmaArray = [sigmaArray; cell2mat(tmpSigma)];
                        tDiff = cell2mat(CombinedResults.IntensityPeaks(indTime,16));
                        tDiff = tDiff(sigNonZInd);


                        xres = (cell2mat(CombinedResults.IntensityPeaks(indTime,14)));
                        yres = (cell2mat(CombinedResults.IntensityPeaks(indTime,15)));
                        tDiffXY = (cell2mat(CombinedResults.IntensityPeaks(indTime,16)));
                        tmpxyres = mean([xres, yres], 2);
                        % remove elements where we have set either x res and y
                        % res to zero due to poor fit.
                        goodInd = find(xres ~= 0 & yres ~= 0);
                        goodInd
                        length(goodInd)
                        length(tmpxyres)
                        tmpxyres = {tmpxyres(goodInd)};
                        tDiffXY = tDiffXY(goodInd);
                        xyres(iT, iR, iField, iPeak, iTime, :) = tmpxyres;
                        xyresArray = [xyresArray; cell2mat(tmpxyres)];
                        % For migration distance, plot vs. EP time.  For sigma
                        % and x-y resolution, plot vs. diffusion time (time in
                        % gel)

                        diffTimeArrayXY =[diffTimeArrayXY; tDiffXY];
                        epTimeArrayMu =[epTimeArrayMu; ...
                            epTimes(iTime)*ones(length(cell2mat(tmpMu)), 1)];
                        diffTimeArraySigma =[diffTimeArraySigma; tDiff];

                        % Old implementation using only EP time (not diffusion
                        % time)
    %                     diffTimeArraySigma =[diffTimeArraySigma; ...
    %                         epTimes(iTime)*ones(length(cell2mat(tmpSigma)), 1)];
    %                     diffTimeArrayXY =[diffTimeArrayXY; ...
    %                         epTimes(iTime)*ones(length(cell2mat(tmpxyres)), 1)];

                        % Potential additional categorization based on
                        % stacks and ROIs for future
    %                     gelReps = unique(cell2mat(CombinedResults.IntensityPeaks(indTime,1)));
    %                     % go through all the gel replicates
    %                     for iGel = 1:length(gelReps)
    %                         indGel = find(cell2mat(CombinedResults.IntensityPeaks(indTime,1)) == gelReps(iGel)); % subset of the indices with this gel replicate number
    %                         stacks = unique(cell2mat(CombinedResults.IntensityPeaks(indGel,6))); 
    %                         % go through all the images from this gel
    %                         for iStack = 1:length(stacks)
    %                             indStack = find(cell2mat(CombinedResults.IntensityPeaks(indGel,6)) == stacks(iStack)); % subset of the indices with this image
    %                             rois = unique(cell2mat(CombinedResults.IntensityPeaks(indStack,7))); 
    %                             for iROI = 1:length(rois)
    %                                 % do nothing for now, as we are just
    %                                 % looking at all the ROIs, images, gels
    %                                 % together.
    %                             end
    %                         end
    %                     end
                    end

                    %**********************************************************
                    % plot and fit the migration distance
                    % figure to hold the boxplot
                    if length(hMu) < iPeak
                        length(hMu)
                        hMu(iPeak) = figure; hold on;
                    else
                        figure(hMu(iPeak)); hold on;
                    end

                    hold on;
                    % plot the boxplot
                    boxplot(muArray, epTimeArrayMu, 'Notch','on', 'Color', cell2mat(colArray(colIdx)));
                    tmp = mu(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    buf = 0.1*min(tmp);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    % putting a legend on a boxplot seems to be a whole disaster
                    if ~suppressLegend
                        bplots = findall(gca,'Tag','Box');
                        for i = 1:length(bplots)
                            bColor = get(bplots(i), 'Color');

                            colourIdx(i) = find(ismember(cell2mat(colArray),bColor, 'rows'));
                        end
                        hLeg = legend(bplots,legendStr(colourIdx), 'Location','southeast');
                    end

                    xlabel('Electrophoresis time (s)');
                    ylabel('Migration distance (\mum)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    title([cell2mat(peaks(iPeak)), ' migration distance'])
                    saveas(hMu(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistance.fig'], 'fig') 
                    saveas(hMu(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistance.svg'], 'svg')
                    saveas(hMu(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistance.png'], 'png')

                    % plot all points with the fit
                    % figure to hold the plot of all points with fit

                    if length(hMuFit) < iPeak
                        length(hMuFit)
                        hMuFit(iPeak) = figure; hold on;
                    else
                        figure(hMuFit(iPeak)); hold on;
                    end

                    CombinedResults.epTimeArrayMu(iT, iR, iField, iPeak) = {epTimeArrayMu};
                    CombinedResults.muArray(iT, iR, iField, iPeak) = {muArray};

                    % plot all the data
                    plot(epTimeArrayMu, muArray, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', 'none', 'Marker', '.');

                    % Add legend entry for the plot 
                    LegendStrFitMu(iPeak, :) = {[LegendStrFitMu{iPeak, :}, {[num2str(gelTs(iT)), ...
                        '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm']}]};

                    if length(epTimes)>1 % only do the fit if we have more than one EP time    
                        % linear fit
                        c = polyfit(epTimeArrayMu,muArray,1);

                        %Add legend entry for the fit
                        LegendStrFitMu(iPeak, :) = {[LegendStrFitMu{iPeak, :}, {[num2str(gelTs(iT)), ...
                            '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm', ...
                            ': y = ', num2str(c(1)), '*x + ', num2str(c(2))]}]};

                        % Plot the fit line
                        y_est = polyval(c,epTimeArrayMu);
                        % Add trend line to plot
                        plot(epTimeArrayMu, y_est, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', ':');
                    end

                    % Add the legend
                    tmp = LegendStrFitMu(iPeak, :);
                    legend(tmp{:});
                    % Save the figure
                    xlabel('Electrophoresis time (s)');
                    ylabel('Migration distance (\mum)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    tbuf = 0.1*min(epTimeArrayMu);
                    if ~isempty(epTimeArrayMu)
                        xlim([min(epTimeArrayMu)-tbuf,max(epTimeArrayMu)+tbuf]);
                    end
                    tmp = mu(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    buf = 0.1*min(tmp);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    title([cell2mat(peaks(iPeak)), ' migration distance'])
                    saveas(hMuFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistanceFit.fig'], 'fig') 
                    saveas(hMuFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistanceFit.svg'], 'svg')
                    saveas(hMuFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_MigrationDistanceFit.png'], 'png')

                    %**********************************************************
                    % plot and fit the peak width
                    if length(hSigma) < iPeak
                        hSigma(iPeak) = figure; hold on;
                    else
                        figure(hSigma(iPeak)); hold on;
                    end
                    hold on;
                    boxplot(sigmaArray.^2, diffTimeArraySigma, 'Notch','on', 'Color', cell2mat(colArray(colIdx)));
                    tmp = sigma(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    tmp = tmp.^2;
                    buf = 0.1*min(tmp);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    % putting a legend on a boxplot seems to be a whole disaster
                    if ~suppressLegend
                        bplots = findall(gca,'Tag','Box');
                        for i = 1:length(bplots)
                            bColor = get(bplots(i), 'Color');

                            colourIdx(i) = find(ismember(cell2mat(colArray),bColor, 'rows'));
                        end
                        hLeg = legend(bplots,legendStr(colourIdx), 'Location','southeast');
                    end

                    xlabel('Electrophoresis time (s)');
                    ylabel('Squared peak width (\mum^2)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    title([cell2mat(peaks(iPeak)), ' peak width'])
                    saveas(hSigma(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidth.fig'], 'fig') 
                    saveas(hSigma(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidth.svg'], 'svg')
                    saveas(hSigma(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidth.png'], 'png')

                    % plot all points with the fit
                    % figure to hold the plot of all points with fit

                    if length(hSigmaFit) < iPeak
                        length(hSigmaFit)
                        hSigmaFit(iPeak) = figure; hold on;
                    else
                        figure(hSigmaFit(iPeak)); hold on;
                    end

                    CombinedResults.sigmaArray(iT, iR, iField, iPeak) = {sigmaArray};

                    % plot all the data
                    plot(diffTimeArraySigma, sigmaArray.^2, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', 'none', 'Marker', '.');

                    % Add legend entry for the plot 
                    LegendStrFitSigma(iPeak, :) = {[LegendStrFitSigma{iPeak, :}, {[num2str(gelTs(iT)), ...
                        '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm']}]};

                    if length(epTimes)>1
                        % linear fit
                        c = polyfit(diffTimeArraySigma,sigmaArray.^2,1);

                        % Add legend entry for the fit
                        LegendStrFitSigma(iPeak, :) = {[LegendStrFitSigma{iPeak, :}, {[num2str(gelTs(iT)), ...
                            '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm', ...
                            ': y = ', num2str(c(1)), '*x + ', num2str(c(2))]}]};

                        % Plot the fit line
                        y_est = polyval(c,diffTimeArraySigma);
                        % Add trend line to plot
                        plot(diffTimeArraySigma, y_est, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', ':');
                    end

                    % Add the legend
                    tmp = LegendStrFitSigma(iPeak, :);
                    legend(tmp{:});
                    % Save the figure
                    xlabel('Diffusion time (s)');
                    ylabel('Squared peak width (\mum^2)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    tbuf = 0.1*min(diffTimeArraySigma);
                    if ~isempty(diffTimeArraySigma)
                        xlim([min(diffTimeArraySigma)-tbuf,max(diffTimeArraySigma)+tbuf]);
                    end
                    tmp = sigma(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    tmp = tmp.^2;
                    buf = 0.1*min(tmp);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    title([cell2mat(peaks(iPeak)), ' peak width'])
                    saveas(hSigmaFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidthFit.fig'], 'fig') 
                    saveas(hSigmaFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidthFit.svg'], 'svg')
                    saveas(hSigmaFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_PeakWidthFit.png'], 'png')

                    %**********************************************************
                    % plot and fit the x-y resolution
                    if length(hXYRes) < iPeak
                        hXYRes(iPeak) = figure; hold on;
                    else
                        figure(hXYRes(iPeak)); hold on;
                    end
                    hold on;
                    boxplot(xyresArray.^2, diffTimeArrayXY, 'Notch','on', 'Color', cell2mat(colArray(colIdx)));
                    tmp = xyres(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    tmp = tmp.^2;
                    buf = 0.1*min(tmp);
                    disp(['min is ', num2str(min(tmp))]);
                    disp(['max is ', num2str(max(tmp))]);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    % putting a legend on a boxplot seems to be a whole disaster
                    if ~suppressLegend
                        bplots = findall(gca,'Tag','Box');
                        for i = 1:length(bplots)
                            bColor = get(bplots(i), 'Color');

                            colourIdx(i) = find(ismember(cell2mat(colArray),bColor, 'rows'));
                        end
                        hLeg = legend(bplots,legendStr(colourIdx), 'Location','southeast');
                    end

                    xlabel('Diffusion time (s)');
                    ylabel('x-y squared peak width (\mum^2)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    title([cell2mat(peaks(iPeak)), ' x-y resolution'])
                    saveas(hXYRes(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_xyRes.fig'], 'fig') 
                    saveas(hXYRes(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_xyRes.svg'], 'svg')
                    saveas(hXYRes(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_xyRes.png'], 'png')

                    % plot all points with the fit
                    % figure to hold the plot of all points with fit

                    if length(hXYResFit) < iPeak
                        length(hXYResFit)
                        hXYResFit(iPeak) = figure; hold on;
                    else
                        figure(hXYResFit(iPeak)); hold on;
                    end

                    CombinedResults.xyresArray(iT, iR, iField, iPeak) = {xyresArray};

                    % plot all the data
                    plot(diffTimeArrayXY, xyresArray.^2, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', 'none', 'Marker', '.');

                    % Add legend entry for the plot
                    LegendStrFitXY(iPeak, :) = {[LegendStrFitXY{iPeak, :}, {[num2str(gelTs(iT)), ...
                        '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm']}]};

                    if length(epTimes)>1
                        % linear fit
                        c = polyfit(diffTimeArrayXY,xyresArray.^2,1);


                        % Add legend entry for the fit
                        LegendStrFitXY(iPeak, :) = {[LegendStrFitXY{iPeak, :}, {[num2str(gelTs(iT)), ...
                            '%T, ',num2str(gelRs(iR)), '%R, ', num2str(eFields(iField)), ' V/cm', ...
                            ': y = ', num2str(c(1)), '*x + ', num2str(c(2))]}]};

                        % Plot the fit line
                        y_est = polyval(c,diffTimeArrayXY);
                        % Add trend line to plot
                        plot(diffTimeArrayXY, y_est, 'Color', cell2mat(colArray(colIdx)), 'LineStyle', ':');
                    end

                    % Add the legend
                    tmp = LegendStrFitXY(iPeak, :);
                    legend(tmp{:});
                    % Save the figure
                    xlabel('Electrophoresis time (s)');
                    ylabel('Squared x-y peak width (\mum^2)');
                    set(gcf, 'Position', [0, 0, figSzX, figSzY])
                    tbuf = 0.1*min(diffTimeArrayXY);
                    if ~isempty(diffTimeArrayXY)
                        xlim([min(diffTimeArrayXY)-tbuf,max(diffTimeArrayXY)+tbuf]);
                    end
                    tmp = xyres(:, :, :, iPeak, :, :);
                    tmp = cell2mat(tmp(:));
                    tmp = tmp.^2;
                    buf = 0.1*min(tmp);
                    if ~isempty(tmp)
                        ylim([min(tmp)-buf,max(tmp)+buf]);
                    end
                    title([cell2mat(peaks(iPeak)), ' x-y resolution'])
                    saveas(hXYResFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_XYResFit.fig'], 'fig') 
                    saveas(hXYResFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_XYResFit.svg'], 'svg')
                    saveas(hXYResFit(iPeak), [StackResults.outputFilenameBase, '_', cell2mat(peaks(iPeak)), '_XYResFit.png'], 'png')

                end
            end
        end
    end
end

    
% plot migration distance, x-y peak width,  and z peak width vs. electrophoresis time as boxplots
% save the images
% average the migration distance, x-y peak width,  and z peak width, in all ROIs of the same image and find the stdev
% average the migration distance, x-y peak width,  and z peak width, in all images using the same electrophoresis conditions.  Use error propagation to calculate the error.
% plot migration distance, x-y peak width,  and z peak width vs. electrophoresis time with error bars
% save the images
% save the data in the master data structure
CombinedResults.Mu = mu;
CombinedResults.Sigma = sigma;
CombinedResults.xyRes = xyres;

end