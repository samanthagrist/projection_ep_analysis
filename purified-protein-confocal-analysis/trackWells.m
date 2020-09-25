%*************************************************************************%
%**          FUNCTION: FIND ROIs BY SEGMENTING EXAMPLE IMAGE            **%
%*************************************************************************%
% Samantha M. Grist
% v.0.1 - 2018-08-03

function StackResults = trackWells(StackResults)
    % initialize the reader and variables
    pauseLength = 0.2; % length of the pause between showing each segmentation
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
    % loop through until we don't find the dark circles anymore,
    % starting at the top of the stack  
    iZ = 0;
    endFlag = 0; % set to 1 when we don't find circles anymore
    startFlag = 0; % set to 1 the first time we find circles (to account for the case where there is stack above the gel top)
    while(~endFlag)
        % increment the slice number 
        iZ = iZ+1;
        
        % read in the image
        A = uint16(zeros(size(I)));
        % if we have defined the example protein, use that.
        if StackResults.exampleProtein
            % calculate the slice number depending on whether we are flipping
            % the stack
            proteinIdx = StackResults.exampleProtein;
            if StackResults.flipStack
                sliceNum = totalZ-iZ + 1;
            else
                sliceNum = iZ;
            end
            disp(['Slice ', num2str(sliceNum)]);
            iPlane = reader.getIndex(sliceNum-1, proteinIdx-1, 0) + 1;
            A = bfGetPlane(reader, iPlane);
        else % sum up the colour channels to create overlay image
            for proteinIdx = 1:length(proteins)
            %for proteinIdx = StackResults.exampleProtein
                % calculate the slice number depending on whether we are flipping
                % the stack
                if StackResults.flipStack
                    sliceNum = totalZ-iZ + 1;
                else
                    sliceNum = iZ;
                end
                disp(['Slice ', num2str(sliceNum)]);
                iPlane = reader.getIndex(sliceNum-1, proteinIdx-1, 0) + 1;
                A = A + bfGetPlane(reader, iPlane);
            end
        end

        % segment the image, looking for dark circles
        StackResults.A = A;
        StackResults.AScaled = sgautoscale(A);
        [StackResults, centroids] = segmentImageDark(StackResults);

        % if we find circles
        if ~isempty(centroids)
            numCentroids(iZ) = length(centroids);
            centroidsArray(iZ) = {centroids};
            % store circle locations in the list of ROI locations
            % (overwrite each loop with the locations from Z plane with most found wells)
            [~, i] = max(numCentroids);
            StackResults.wellLocs = centroidsArray(i);
            % set the startFlag so we know we have found wells
            startFlag = 1; 
        else
            numCentroids(iZ) = 0;
            %centroidsArray
            centroidsArray(iZ) = {0};
            if startFlag % if we don't find circles and startFlag == 1
            % set endFlag to 1
            endFlag = 1;
            end
            % if we don't find circles but startFlag == 0, do nothing
        end 
        pause(pauseLength);
    end

    % the ending iZ value is the first one in which we no longer find
    % dark circles
    StackResults.topZ = iZ;
    disp(['topZ = ', num2str(iZ)]);
    disp(['found ', num2str(numCentroids(i)), ' wells']);
end