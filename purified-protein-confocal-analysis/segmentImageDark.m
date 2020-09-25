%*************************************************************************%
%**      FUNCTION: FIND ROIs BY SEGMENTING EXAMPLE IMAGE FOR WELLS      **%
%*************************************************************************%
% Samantha M. Grist
% v.0.1 - 2018-08-03

function [StackResults, centroids] = segmentImageDark(StackResults)
    % tuning parameters
    gaussfiltSz = 3; % sigma for Gaussian blur filter
    threshSensitivity = 0.81; %0.79; % sensitivity for image adaptive thresholding
    closeDiskSz = 3; % radius for disk used in imclose operation (filling holes in segmented wells)

    % median filter to improve segmentation results
    A = medfilt2(StackResults.A);
    AScaled = medfilt2(StackResults.AScaled);
   
    % gaussian filter to blur the noise
    Agauss = imgaussfilt(A, gaussfiltSz);
    
    % threshold the image using adaptive thresholding
    BWtest = imbinarize(Agauss, 'adaptive', 'Sensitivity', threshSensitivity);
    
    % invert the thresholded image (making the wells bright)
    BWC = ~BWtest;
    
    % morphologic close to fill in the wells 
    se = strel('disk', closeDiskSz); BWCclose = imclose(BWC, se);
    
    % remove objects smaller than circle area corresponding to 0.7 time the estimated radius
    estArea = double(uint16((0.7*StackResults.estRadius)^2*pi));
    BWCfinal = bwareaopen(BWCclose, estArea);
    
    % find dark circles in image
    s = regionprops(BWCfinal,'centroid');
    centroids = cat(1, s.Centroid);
    
    figure(1); imshow(3*AScaled);
    if ~isempty(centroids)
        figure(1); hold on;
        plot(centroids(:,1),centroids(:,2), 'r*')
    end
end
