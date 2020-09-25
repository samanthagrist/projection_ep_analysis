%**                                                                     **%
%**                       ANALYZE LYSIS DATA                            **%
%**                        October 9th, 2018                            **%
%**                        Samantha M. Grist                            **%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%
close all;
clear;

% set the root directory where the image analysis scripts and functions live, set up so that it runs on a Windows and MacOS system with different paths
if ispc % Windows path
    dataAnalysisRoot = 'path1';
    trackingRoot = [dataAnalysisRoot,'\tracking'];
else % MacOS path
    dataAnalysisRoot = 'path2';
    trackingRoot = [dataAnalysisRoot,'/tracking'];
end

addpath(dataAnalysisRoot);
addpath(trackingRoot);

configFilename = 'Initialization.xlsx';
dataFilename = 'output.mat';
restart = 0; %if set to 1, do not load previous data (will save over if filename is the same)
dbg = 0;
figNum = 1;

if exist(dataFilename, 'file') && ~restart
    load(dataFilename);
    % Read in the stage location configuration information from a .csv file   
    ImgInfoTable = readtable(configFilename);
    ImgInfo = table2cell(ImgInfoTable);
    maxNumFiles = size(ImgInfo,1); % the number of files in the current directory gives an upper end of the images we can analyze.
    maxMonitoringLength = 300; % max length (in pixels) of the monitoring stacks
    restart = 0; %reset the restart flag
else
    % Read in the stage location configuration information from a .csv file
    ImgInfoTable = readtable(configFilename);
    ImgInfo = table2cell(ImgInfoTable);
    maxNumFiles = size(ImgInfo,1); % the number of files in the current directory gives an upper end of the images we can analyze.
    maxMonitoringLength = 300; % max length (in pixels) of the monitoring stacks
    analyzedFlag = zeros(maxNumFiles, 1);
end

if restart
    analyzedFlag = zeros(maxNumFiles, 1);
end
% loop through all the files, analyzing the lysis data
for fileIdx = 1:maxNumFiles
    disp(['*****Now analyzing file ', num2str(fileIdx)]);
    conditionNum(fileIdx) = ImgInfo{fileIdx, 1};
    filename = ImgInfo{fileIdx, 2};
    outputFilenameRoot = ImgInfo{fileIdx, 3};
    outputFilename = [outputFilenameRoot, '.avi'];
    conditionText(fileIdx) = ImgInfo(fileIdx, 4);
    startIdx = ImgInfo{fileIdx, 5};
    fudge = ImgInfo{fileIdx, 11};
    if startIdx > 1
        excludeFrames = [1:startIdx];
    else
        excludeFrames = [];
    end
    tInterval(fileIdx) = ImgInfo{fileIdx, 6};
    pixScale(fileIdx) = ImgInfo{fileIdx, 7};
    repNum(fileIdx) = ImgInfo{fileIdx, 8};
    integrationRadius(fileIdx) = ImgInfo{fileIdx, 9};
    wellRadius(fileIdx) = ImgInfo{fileIdx, 10};
    figNum = 1;
    
    if ~analyzedFlag(fileIdx)
        [tArray{conditionNum(fileIdx), repNum(fileIdx)}, normParticleIntensity{conditionNum(fileIdx), repNum(fileIdx)}, ...
            wellIntensity{conditionNum(fileIdx), repNum(fileIdx)}, maxIntensity{conditionNum(fileIdx), repNum(fileIdx)}, ...
            normWellIntensity{conditionNum(fileIdx), repNum(fileIdx)}, normMaxIntensity{conditionNum(fileIdx), repNum(fileIdx)}, ...
            proteinWidth{conditionNum(fileIdx), repNum(fileIdx)}] = ...
            AnalyzeLysisData(filename, outputFilename, excludeFrames, ...
            tInterval(fileIdx), integrationRadius(fileIdx), wellRadius(fileIdx), dbg, fudge, pixScale(fileIdx), figNum);
        analyzedFlag(fileIdx) = 1;
        save(dataFilename);
    end

    conditionTxtRearr(conditionNum(fileIdx), repNum(fileIdx)) = conditionText(fileIdx);
end
save(dataFilename);

%% plot the datas
cols = distinguishable_colors(max(conditionNum)+1);
% plot the analyzed lysis data
for conditionIdx = 1:max(conditionNum)
    hParticle = figure(5+conditionIdx); hold on;
    hWell = figure(5+conditionIdx+max(conditionNum)); hold on;
    hMax = figure(5+conditionIdx+2*max(conditionNum)); hold on;
    hWidth = figure(5+conditionIdx+3*max(conditionNum)); hold on;
    for repIdx = 1:max(repNum)
        if repIdx == 1
            disp(['ConditionIdx = ' , num2str(conditionIdx), ', repIdx = ', num2str(repIdx)]);
            % plot the lines on each plot, making the first legend entry visible
            
            % store the four types of data in temporary variables
            tempParticle = normParticleIntensity{conditionIdx, repIdx};
            tempWell = normWellIntensity{conditionIdx, repIdx};
            tempMax = normMaxIntensity{conditionIdx, repIdx};
            tempWidth = proteinWidth{conditionIdx, repIdx};
            % loop through all the tracked particles
            particleTracker = 0;
            for particleIdx = 1:size(tempParticle, 2)
                if any(tempParticle(:,particleIdx))
                    particleTracker = particleTracker + 1;
                    if particleTracker == 1  % plot the first line as visible legend entry
                        disp(['First particle - adding to legend entry'])
                        tempTime = tArray{conditionIdx, repIdx};
                        tempTime = tempTime - tempTime(1);
                        figure(2); hold on;
                        plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                            nonzeros(tempParticle(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(3); hold on;
                        plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                            nonzeros(tempWell(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(4); hold on;
                        plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                            nonzeros(tempMax(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(5); hold on;
                        plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                            nonzeros(tempWidth(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        % Plot the individual figures too 
                        figure(hParticle); hold on;
                        plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                            nonzeros(tempParticle(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(hWell); hold on;
                        plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                            nonzeros(tempWell(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(hMax); hold on;
                        plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                            nonzeros(tempMax(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        figure(hWidth); hold on;
                        plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                            nonzeros(tempWidth(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','on');
                        % store the legend entry
                        legendStr(conditionIdx) = {strrep(conditionTxtRearr{conditionIdx, 1}, '_', ' ')};
                    else  % make all subsequent lines not visible in legend
                        tempTime = tArray{conditionIdx, repIdx};
                        tempTime = tempTime - tempTime(1);
                        figure(2); hold on;
                        plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                            nonzeros(tempParticle(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(3); hold on;
                        plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                            nonzeros(tempWell(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(4); hold on;
                        plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                            nonzeros(tempMax(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(5); hold on;
                        plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                            nonzeros(tempWidth(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                         
                         % Plot the individual figures too 
                        figure(hParticle); hold on;
                        plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                            nonzeros(tempParticle(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(hWell); hold on;
                        plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                            nonzeros(tempWell(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(hMax); hold on;
                        plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                            nonzeros(tempMax(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                        figure(hWidth); hold on;
                        plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                            nonzeros(tempWidth(:, particleIdx)), ...
                             'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                    end
                end
            end
            
            
        else
            % plot the lines on each plot, making the legend entry invisible
            % store the four types of data in temporary variables
            tempParticle = normParticleIntensity{conditionIdx, repIdx};
            tempWell = normWellIntensity{conditionIdx, repIdx};
            tempMax = normMaxIntensity{conditionIdx, repIdx};
            tempWidth = proteinWidth{conditionIdx, repIdx};
            % loop through all the tracked particles
            for particleIdx = 1:min([size(tempParticle, 2), size(tempWell, 2), size(tempMax, 2), size(tempWidth, 2)])
                disp(['ParticleIdx = ', num2str(particleIdx)]);
                tempTime = tArray{conditionIdx, repIdx};
                tempTime = tempTime - tempTime(1);
                figure(2); hold on;
                plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                    nonzeros(tempParticle(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(3); hold on;
                plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                    nonzeros(tempWell(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(4); hold on;
                plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                    nonzeros(tempMax(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(5); hold on;
                particleIdx
                size(tempTime)
                size(tempWidth)
                plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                    nonzeros(tempWidth(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                 
                 % Plot the individual figures too 
                figure(hParticle); hold on;
                plot(tempTime(1:length(nonzeros(tempParticle(:, particleIdx)))), ...
                    nonzeros(tempParticle(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(hWell); hold on;
                plot(tempTime(1:length(nonzeros(tempWell(:, particleIdx)))), ...
                    nonzeros(tempWell(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(hMax); hold on;
                plot(tempTime(1:length(nonzeros(tempMax(:, particleIdx)))), ...
                    nonzeros(tempMax(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
                figure(hWidth); hold on;
                plot(tempTime(1:length(nonzeros(tempWidth(:, particleIdx)))), ...
                    nonzeros(tempWidth(:, particleIdx)), ...
                     'Color', cols(conditionIdx, :), 'HandleVisibility','off');
            end
        end  
    end
    
    % add the legends and axis labels and save each individual condition figure
    h = figure(hParticle); title(legendStr(conditionIdx)); xlabel('time (s)'); ylabel('Normalized integrated particle intensity');
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normParticleIntensity.fig']);
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normParticleIntensity.svg']);
    h = figure(hWell); title(legendStr(conditionIdx)); xlabel('time (s)'); ylabel('Normalized integrated particle intensity in well');
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normWellIntensity.fig']);
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normWellIntensity.svg']);
    h = figure(hMax); title(legendStr(conditionIdx)); xlabel('time (s)'); ylabel('Normalized maximum particle intensity');
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normMaxIntensity.fig']);
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_normMaxIntensity.svg']);
    h = figure(hWidth); title(legendStr(conditionIdx)); xlabel('time (s)'); ylabel('Protein peak width \sigma (\mum)');
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_peakWidth.fig']);
    saveas(h, [strrep(dataFilename, '.mat', ''), '_condition', num2str(conditionIdx), '_peakWidth.svg']);
end

% add the legends and axis labels and save each overlay figure
h = figure(2); legend(legendStr); xlabel('time (s)'); ylabel('Normalized integrated particle intensity');
saveas(h, [strrep(dataFilename, '.mat', ''), '_normParticleIntensity.fig']);
saveas(h, [strrep(dataFilename, '.mat', ''), '_normParticleIntensity.svg']);
h = figure(3); legend(legendStr); xlabel('time (s)'); ylabel('Normalized integrated particle intensity in well');
saveas(h, [strrep(dataFilename, '.mat', ''), '_normWellIntensity.fig']);
saveas(h, [strrep(dataFilename, '.mat', ''), '_normWellIntensity.svg']);
h = figure(4); legend(legendStr); xlabel('time (s)'); ylabel('Normalized maximum particle intensity');
saveas(h, [strrep(dataFilename, '.mat', ''), '_normMaxIntensity.fig']);
saveas(h, [strrep(dataFilename, '.mat', ''), '_normMaxIntensity.svg']);
h = figure(5); legend(legendStr); xlabel('time (s)'); ylabel('Protein peak width \sigma (\mum)');
saveas(h, [strrep(dataFilename, '.mat', ''), '_peakWidth.fig']);
saveas(h, [strrep(dataFilename, '.mat', ''), '_peakWidth.svg']);

