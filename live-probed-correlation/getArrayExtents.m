
function StackResults = getArrayExtents(StackResults, stackIdx, ADisp, figNum)
    wLine = 4;
    dotSz = 3;
    % Display image
    h = figure(figNum); imshow(ADisp);
    
    % Store the top left well position
    title('Zoom into the bottom right corner then press any key')
    zoom on
    while ~waitforbuttonpress()
    end
    title('Click in the bottom right corner of the gel (opposite the cut corner)')
    [x,y,button] = ginput(1)
    StackResults(stackIdx).arrayPos{figNum,1} = [x, y];
    % Burn the top left position into the display image
    ADisp = insertShape(ADisp, 'FilledCircle', [x,y,dotSz], 'Color', 'r');
    pos = get(h, 'Position');
    figure(h); imshow(ADisp), 
    set(h, 'Position', pos);
    
    % Store the top right well position
    title('Zoom into the bottom left corner then press any key')
    zoom on
    while ~waitforbuttonpress()
    end
    title('Click in the bottom left corner of the gel')
    [x,y,button] = ginput(1)
    StackResults(stackIdx).arrayPos{figNum,2} = [x, y];
    % Burn the top right position into the display image
    ADisp = insertShape(ADisp, 'FilledCircle', [x,y,dotSz], 'Color', 'r');
    pos = get(h, 'Position');
    figure(h); imshow(ADisp), 
    set(h, 'Position', pos);
    
    % Store the bottom left well position
    title('Zoom into the top right corner then press any key')
    zoom on
    while ~waitforbuttonpress()
    end
    title('Click in the top right corner of the gel')
    [x,y,button] = ginput(1)
    StackResults(stackIdx).arrayPos{figNum, 3} = [x, y];
    % Burn the bottom left position into the display image
    ADisp = insertShape(ADisp, 'FilledCircle', [x,y,dotSz], 'Color', 'r');
    % Burn lines showing the clicked positions into the image and replot
    ADisp = insertShape(ADisp, 'Line', ...
        [StackResults(stackIdx).arrayPos{figNum,1},StackResults(stackIdx).arrayPos{figNum,2}], ...
        'LineWidth', wLine, 'Color', 'w');
    ADisp = insertShape(ADisp, 'Line', ...
        [StackResults(stackIdx).arrayPos{figNum,1},StackResults(stackIdx).arrayPos{figNum,3}], ...
        'LineWidth', wLine, 'Color', 'w');
    pos = get(h, 'Position');
    figure(h); imshow(ADisp), 
    set(h, 'Position', pos);
    
    % log that we have found the extents
    StackResults(stackIdx).extentsFound(figNum) = 1;
end
