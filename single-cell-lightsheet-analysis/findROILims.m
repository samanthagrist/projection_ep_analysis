% find the start and end indices for the ROI of interest, and check that
% they aren't too close to the image edge
% v0.1 Samantha Grist 2019-08-07
function [xStart, yStart, xEnd, yEnd] = findROILims(StackResults, roiIdx, roiSz, ASz, suppressOut)
                    
    yStart = uint16(StackResults.ROILocX(roiIdx));
    yEnd = uint16(yStart+roiSz);
    xStart = uint16(StackResults.ROILocY(roiIdx));
    xEnd = uint16(xStart+roiSz);
    % check the values
    if xStart < 1
        disp(['[findROILims] Before correction xStart: ', num2str(xStart), ', xEnd: ', num2str(xEnd), ', yStart: ', num2str(yStart), ', yEnd: ', num2str(yEnd)])
        disp(['[findROILims] Adjusting ROI ', num2str(roiIdx), ' xStart']);
        xEnd = xEnd + (1-xStart);
        xStart = xStart + (1-xStart);
    end
    if yStart < 1
        disp(['[findROILims] Before correction xStart: ', num2str(xStart), ', xEnd: ', num2str(xEnd), ', yStart: ', num2str(yStart), ', yEnd: ', num2str(yEnd)])
        disp(['[findROILims] Adjusting ROI ', num2str(roiIdx), ' yStart']);
        yEnd = yEnd + (1-yStart);
        yStart = yStart + (1-yStart);
    end
    if xEnd > ASz(1)
        disp(['[findROILims] Before correction xStart: ', num2str(xStart), ', xEnd: ', num2str(xEnd), ', yStart: ', num2str(yStart), ', yEnd: ', num2str(yEnd)])
        disp(['[findROILims] Adjusting ROI ', num2str(roiIdx), ' xEnd']);
        xStart = xStart - (xEnd - ASz(1));
        xEnd = xEnd - (xEnd - ASz(1));
    end
    if yEnd > ASz(2)
        disp(['[findROILims] Before correction xStart: ', num2str(xStart), ', xEnd: ', num2str(xEnd), ', yStart: ', num2str(yStart), ', yEnd: ', num2str(yEnd)])
        disp(['[findROILims] Adjusting ROI ', num2str(roiIdx), ' yEnd']);
        yStart = yStart - (yEnd - ASz(2));
        yEnd = yEnd - (yEnd - ASz(2));
    end
    if ~suppressOut
        disp(['[findROILims] xStart: ', num2str(xStart), ', xEnd: ', num2str(xEnd), ', yStart: ', num2str(yStart), ', yEnd: ', num2str(yEnd)])
    end

end