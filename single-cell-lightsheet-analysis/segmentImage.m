%*************************************************************************%
%**          FUNCTION: FIND ROIs BY SEGMENTING EXAMPLE IMAGE            **%
%*************************************************************************%
% Samantha M. Grist
% v.0.3 - 2019-08-06

function StackResults = segmentImage(StackResults)
    dbg = 0;
    
    A = StackResults.ASum;
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
    A = imgaussfilt(A, 5);
   
    % Adaptive thresholding
    thresh = graythresh(A);
    %StackResults.thresholdFudge
    ABW = imbinarize(A, 'adaptive', 'sensitivity', StackResults.thresholdFudge);
    %ABW = im2bw(A, StackResults.thresholdFudge*thresh);
    if dbg
        figure; imshow(ABW); pause(1);
    end
    se = strel('disk', 5);
    ABW = imclose(ABW, se);
    se = strel('disk', 7);
    ABW = imopen(ABW, se);
    if dbg
        figure; imshow(ABW); pause(1);
    end
    ABW = imfill(ABW, 'holes');
    ABW = imclose(ABW, se);
    se = strel('disk', 10);
    ABW = imopen(ABW, se);
    ABW = imclose(ABW, se);
    if dbg
        figure; imshow(ABW); pause(1);
    end
    se = strel('disk', 15);
    ABW = imopen(ABW, se);
    ABW=imclearborder(ABW);
    figure(2); imshow(ABW);
    
    figure(1); imshow(sgautoscale(A)); 
    
    
    s = regionprops(ABW,'centroid');
    centroids = cat(1, s.Centroid);
    if length(centroids)      
        % use the centroids to compute and draw the ROIs
        xmin = centroids(:,1)-StackResults.roiSz/2;
        ymin = centroids(:,2)-StackResults.roiSz/2;
        StackResults.centroids = {centroids};
        StackResults.ROILocX = xmin;
        StackResults.ROILocY = ymin;
        AROIs = sgautoscale(A);
        StackResults.numROIs = length(xmin);
        for roiIdx = 1:length(xmin)
            AROIs = insertShape(AROIs, 'Rectangle', [xmin(roiIdx), ymin(roiIdx), StackResults.roiSz, StackResults.roiSz], 'Color', 'r', 'Linewidth', 2);
        end
        % plot the ROIs and segmented centroids
        figure(250); imshow(AROIs); hold on;
        plot(centroids(:,1),centroids(:,2), 'r*')
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

