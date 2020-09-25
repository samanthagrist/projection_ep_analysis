%*************************************************************************%
%**          FUNCTION: FIND ROIs BY SEGMENTING EXAMPLE IMAGE            **%
%*************************************************************************%
% Samantha M. Grist
% v.0.3 - 2019-08-06

function StackResults = segmentImage(StackResults, A, sens, minSize, edgeBuffer)
    dbg = 0;
    
    % if an RGB image was passed in, change it to grayscale by averaging in
    % the third dimension
    if length(size(A)) > 2
        if length(size(A)) > 3
            error('Image matrix has too many dimensions - only 3 dimensions are supported')
        else
            Atemp = double(mean(A, 3));
            A = Atemp./max(max(Atemp));
        end
    else
        A = A./max(max(A));
    end
    %A = medfilt2();
    A = imgaussfilt(A, 10);
   
    % Adaptive thresholding
    T = adaptthresh(A, sens, 'NeighborhoodSize', [2*StackResults.roiSz+1,StackResults.roiSz+1]);
    ABW = imbinarize(A, T);
    if dbg
        figure; imshow(ABW); pause(1);
    end
    % remove objects near the border
    ASz = size(ABW);
    ABW = ABW(edgeBuffer:ASz(1)-edgeBuffer, edgeBuffer:ASz(2)-edgeBuffer);
    
    
    if dbg
        figure; imshow(ABW); pause(1);
    end
    se = strel('disk', 2);
    ABW = imclose(ABW, se);
    se = strel('disk', 4);
    ABW = imopen(ABW, se);
    if dbg
        figure; imshow(ABW); pause(1);
    end
    ABW=bwareaopen(ABW, minSize);
    ABW=imclearborder(ABW);
    
    ABW2 = zeros(ASz);
    ABW2(edgeBuffer:ASz(1)-edgeBuffer, edgeBuffer:ASz(2)-edgeBuffer) = ABW;
    
    figure(300); imshow(ABW);
    
    figure(301); imshow(sgautoscale(A)); 
    
    
    s = regionprops(ABW,'centroid');
    centroids = cat(1, s.Centroid)+edgeBuffer;
    if length(centroids)      
        % use the centroids to compute and draw the ROIs
        xmin = centroids(:,1);
        ymin = centroids(:,2);
        StackResults.centroids = {centroids};
        StackResults.ROILocX = xmin;
        StackResults.ROILocY = ymin;
        AROIs = sgautoscale(A);
        StackResults.numROIs = length(xmin);
%         for roiIdx = 1:length(xmin)
%             AROIs = insertShape(AROIs, 'Rectangle', [xmin(roiIdx), ymin(roiIdx), StackResults.roiSz, StackResults.roiSz], 'Color', 'r', 'Linewidth', 2);
%         end
        % plot the ROIs and segmented centroids
        figure(250); imshow(AROIs); hold on;
        plot(centroids(:,1),centroids(:,2), '.r')
        disp(['[segmentImage]: ', num2str(StackResults.numROIs), ' ROIs found.']);
        
        pause(0.01);
    
    else
        StackResults.centroids = {[]};
        StackResults.ROILocX = [];
        StackResults.ROILocY = [];
        disp('[segmentImage]: No ROIs found.');
    end
    hold off
end

