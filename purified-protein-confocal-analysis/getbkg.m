%*************************************************************************%
%**        Get image background average intensity and std dev           **%
%*************************************************************************%
% Samantha M. Grist
% v.0.2 - 2018-03-01

function StackResults = getbkg(StackResults)
% Process the image (at the first Z value of TopZ) to find the background 
% region in the gel and calculate the background average and noise for 
% each colour channel

% set up the bioFormats reader and get the number of planes in the stack
% and number of Z slices
reader = StackResults.imgReader;
reader.setSeries(0);
planes = reader.getImageCount;
numZ = planes/length(StackResults.proteins);
planes


% calculate the z slice we want to analyze for background regions.  Choose
% 2 slices down (into the gel) from the topZ (microwell bottom)
if ~StackResults.flipStack
    zVal = StackResults.topZ + 3;
else
    zVal = numZ-(StackResults.topZ + 3)+1;
end

% loop through the colour channels and sum to find the background regions
for colourIdx = 1:length(StackResults.proteins)
    % set up the plane index from the Z val and colour
    iPlane = reader.getIndex(zVal-1, colourIdx-1, 0) + 1;
    
    % read in the image
    I = bfGetPlane(reader, iPlane);
    
    % add to the summed colour channel image
    if colourIdx == 1
        ISum = I;
    else
        ISum = ISum + I;
    end
end

% intensity threshold segmentation and dilation to find the background
% region
figure(1); imshow(histeq(ISum)); title('Image');
ISum = medfilt2(ISum);
ISum = imgaussfilt(ISum, 5);
thresh = graythresh(ISum);
BW2 = im2bw(ISum, thresh*1.1); 
%BW2 = im2bw(ISum, thresh*StackResults.thresholdFudge); 
figure(3); imshow(BW2); title('B&W image');
se = strel('disk', 2);
BW = imclose(BW2, se);
se = strel('disk', 3);
BW = imopen(BW, se);
BW = imfill(BW, 'holes');
se = strel('disk', 10);
BW = imopen(BW, se);
se = strel('disk', StackResults.bkgDilation);
IBkg = ~imdilate(BW, se);
figure(2); imshow(IBkg); title('Background region'); pause(0.1);
    
% loop through the colour channels to find the background regions
for colourIdx = 1:length(StackResults.proteins)
    % set up the plane index from the Z val and colour
    iPlane = reader.getIndex(zVal-1, colourIdx-1, 0) + 1;
    
    % read in the image
    I = bfGetPlane(reader, iPlane);
    
    % calculate the background mean and standard deviation
    mean2(I(IBkg))
    StackResults.bkgMean(colourIdx) = mean2(I(IBkg));
    StackResults.bkgStd(colourIdx) = std2(I(IBkg));
end
end