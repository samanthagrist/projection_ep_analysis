%*************************************************************************%
%**          FUNCTION: FIND ROIs BY SEGMENTING EXAMPLE IMAGE            **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2018-03-01

function StackResults = segmentImage(StackResults)
    
    A = StackResults.A;
    %A = medfilt2();
    A = imgaussfilt(A, 5);
    AScaled = medfilt2(StackResults.AScaled);
    AScaled = imgaussfilt(AScaled, 5);
   
    % Intensity threshold using Otsu's method to find the threshold limit
    thresh = graythresh(A);
    %StackResults.thresholdFudge
    ABW = imbinarize(A, 'adaptive', 'sensitivity', StackResults.thresholdFudge);
    %ABW = im2bw(A, StackResults.thresholdFudge*thresh);
    %figure; imshow(ABW); pause(1);
    se = strel('disk', 5);
    ABW = imclose(ABW, se);
    se = strel('disk', 7);
    ABW = imopen(ABW, se);
    %figure; imshow(ABW); pause(1);
    ABW = imfill(ABW, 'holes');
    ABW = imclose(ABW, se);
    se = strel('disk', 20);
    ABW = imopen(ABW, se);
    ABW = imclose(ABW, se);
    %figure; imshow(ABW); pause(1);
    se = strel('disk', 25);
    ABW = imopen(ABW, se);
    ABW=imclearborder(ABW);
    figure(2); imshow(ABW);
    
    figure(1); imshow(AScaled); 
    
    
    s = regionprops(ABW,'centroid');
    centroids = cat(1, s.Centroid);
    if length(centroids)
        figure(1); hold on;
        plot(centroids(:,1),centroids(:,2), 'r*')
        pause(0.01);
        StackResults.centroids(StackResults.currProtein,  StackResults.currZ, :) = {centroids};
    
    else
        StackResults.centroids(StackResults.currProtein,  StackResults.currZ, :) = {[]};
    end
    hold off
end
