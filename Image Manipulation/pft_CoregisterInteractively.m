function [ Registered, RowShift, ColShift, PlnShift ] = pft_CoregisterInteractively(Moving, Fixed, DR, DC, DP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fetch some image dimensions first
DIMS = size(Fixed);

ROWS = DIMS(1); 
COLS = DIMS(2); 
PLNS = DIMS(3); 

% Size the display axes
Rows = 3*DIMS(1); 
Cols = 3*DIMS(2); 
Plns = 3*DIMS(3); 

% Initialize the slider values - take care of odd/even dimensions with one instruction
Row = floor((ROWS + 1)/2);  
Col = floor((COLS + 1)/2);
Pln = floor((PLNS + 1)/2);

RowShift = DR;
ColShift = DC;
PlnShift = DP;

% Provide a zero-padding buffer for later circular-shifted images
Pads = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Acquire some screen dimensions
area = get(0, 'ScreenSize');

% This is the position for the GUI figure as a whole
rect = zeros([1, 4], 'double');

rect(3) = area(3) - 256;
rect(4) = area(4) - 128;
rect(1) = (area(3) - rect(3))/2;
rect(2) = (area(4) - rect(4))/2;

hf = figure('Name', 'Manual Co-Registration', 'MenuBar', 'none', 'NumberTitle', 'off', 'Renderer', 'OpenGL', 'Position', rect);

set(gcf, 'CloseRequestFcn', @my_closefcn);  % Trap the close button
set(gcf, 'KeyPressFcn', @my_keypressfcn);   % Trap a RETURN or an ESCAPE from the keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position 6 sliders - 3 for the fixed position and 3 for the Moving offset - with a uniform layout
DX = 32;
DY = 4;
WD = rect(3)/2 - 2*DX;
HT = 32;

here = [ DX, 4*DY, WD, HT ];
hPlnSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', Pln, 'Min', Pln-24, 'Max', Pln+24, 'SliderStep', [1.0, 4.0]/48.0, ...
                       'Position', here, 'Callback', @pln_slider_callback);                   
hPlnSliderListener = addlistener(hPlnSlider, 'ContinuousValueChange', @continuous_pln_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hPlnSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Slice: %1d', Pln), 'FontSize', 16, 'FontWeight', 'bold');

here = here + [ 0, HT + DY, 0, 0 ];
hColSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', Col, 'Min', Col-48, 'Max', Col+48, 'SliderStep', [1.0, 4.0]/96.0, ...
                       'Position', here, 'Callback', @col_slider_callback);                   
hColSliderListener = addlistener(hColSlider, 'ContinuousValueChange', @continuous_col_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hColSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Column: %1d', Pln), 'FontSize', 16, 'FontWeight', 'bold');

here = here + [ 0, HT + DY, 0, 0 ];
hRowSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', Row, 'Min', Row-48, 'Max', Row+48, 'SliderStep', [1.0, 4.0]/96.0, ...
                       'Position', here, 'Callback', @row_slider_callback);                   
hRowSliderListener = addlistener(hRowSlider, 'ContinuousValueChange', @continuous_row_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hRowSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Row: %1d', Row), 'FontSize', 16, 'FontWeight', 'bold');

here = [ DX, 4*DY, WD, HT ];
here = here + [ WD + 2*DX, 0, 0, 0 ];
hPlnShiftSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', PlnShift, 'Min', -24, 'Max', 24, 'SliderStep', [1.0, 4.0]/48.0, ...
                            'Position', here, 'Callback', @pln_shift_slider_callback);                   
hPlnShiftSliderListener = addlistener(hPlnShiftSlider, 'ContinuousValueChange', @continuous_pln_shift_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hPlnShiftSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Slice shift: %1d', PlnShift), 'FontSize', 16, 'FontWeight', 'bold');

here = here + [ 0, HT + DY, 0, 0 ];
hColShiftSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', ColShift, 'Min', -48, 'Max', 48, 'SliderStep', [1.0, 4.0]/96.0, ...
                            'Position', here, 'Callback', @col_shift_slider_callback);                   
hColShiftSliderListener = addlistener(hColShiftSlider, 'ContinuousValueChange', @continuous_col_shift_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hColShiftSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Column shift: %1d', ColShift), 'FontSize', 16, 'FontWeight', 'bold');

here = here + [ 0, HT + DY, 0, 0 ];
hRowShiftSlider = uicontrol('Parent', hf, 'Style', 'Slider', 'Value', RowShift, 'Min', -48, 'Max', 48, 'SliderStep', [1.0, 4.0]/96.0, ...
                            'Position', here, 'Callback', @row_shift_slider_callback);                   
hRowShiftSliderListener = addlistener(hRowShiftSlider, 'ContinuousValueChange', @continuous_row_shift_slider_callback);

here = here + [ 0, HT + DY, 0, 0 ];
hRowShiftSliderText = uicontrol('Parent', hf, 'Style', 'Text', 'Position', here, 'String', sprintf('Row shift: %1d', RowShift), 'FontSize', 16, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position the push-button to commit changes and exit
here = [ DX, 4*DY, WD, HT ];
here = here + [ 0, 8*DY + 6*HT, 2*DX + WD, HT ];
hApplyButton = uicontrol('Parent', hf, 'Style', 'Pushbutton', 'Position', here, 'String', 'Apply Shift and Exit', 'FontSize', 18, 'FontWeight', 'bold', 'Callback', @pushbutton_callback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now position the axes to create 3 orthogonal views of the overlaid i/p arrays
DQ = 4;

here = [ rect(3)/2 - Cols/2, 0, Cols, Rows ];
here = here + [ 0, 8*DY + 9*HT, 0, 0 ];
hAxes12 = axes('Parent', hf, 'Units', 'Pixels', 'Position', here);
imshow(repmat(uint8(0), [Rows, Cols]), [], 'Parent', hAxes12);

here = [ rect(3)/2 - Cols/2, 0, Cols, Rows ];
here = here + [ 0, 8*DY + 9*HT, 0, 0 ];
here(1) = here(1) + DQ + Cols;
here(3) = Plns;
here(4) = Rows;
hAxes13 = axes('Parent', hf, 'Units', 'Pixels', 'Position', here);
imshow(repmat(uint8(64), [Rows, Plns]), [], 'Parent', hAxes13);

here = [ rect(3)/2 - Cols/2, 0, Cols, Rows ];
here = here + [ 0, 8*DY + 9*HT, 0, 0 ];
here(2) = here(2) + DQ + Rows;
here(3) = Cols;
here(4) = Plns;
hAxes23 = axes('Parent', hf, 'Units', 'Pixels', 'Position', here);
imshow(repmat(uint8(128), [Plns, Cols]), [], 'Parent', hAxes23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now initialize the display and block user interaction from the command line while the figure continues to exist (as a 'modal' object)
UpdateDisplay;

uiwait(hf);

    function pln_slider_callback(hObj, eventdata)
        Pln = get(hObj, 'Value');
        Pln = round(Pln);
        str = sprintf('Slice: %1d', Pln);
        set(hPlnSliderText, 'String', str);
        
        UpdateDisplay;         
    end

    function continuous_pln_slider_callback(hObj, eventdata)
        Pln = get(hObj, 'Value');
        Pln = round(Pln);
        str = sprintf('Slice: %1d', Pln);
        set(hPlnSliderText, 'String', str);
        
        UpdateDisplay; 
    end

    function col_slider_callback(hObj, eventdata)
        Col = get(hObj, 'Value');
        Col = round(Col);
        str = sprintf('Column: %1d', Col);
        set(hColSliderText, 'String', str);
        
        UpdateDisplay;         
    end

    function continuous_col_slider_callback(hObj, eventdata)
        Col = get(hObj, 'Value');
        Col = round(Col);
        str = sprintf('Column: %1d', Col);
        set(hColSliderText, 'String', str);
        
        UpdateDisplay;      
    end

    function row_slider_callback(hObj, eventdata)
        Row = get(hObj, 'Value');
        Row = round(Row);
        str = sprintf('Row: %1d', Row);
        set(hRowSliderText, 'String', str);
        
        UpdateDisplay;      
    end

    function continuous_row_slider_callback(hObj, eventdata)
        Row = get(hObj, 'Value');
        Row = round(Row);
        str = sprintf('Row: %1d', Row);
        set(hRowSliderText, 'String', str);
        
        UpdateDisplay;   
    end

    function pln_shift_slider_callback(hObj, eventdata)
        PlnShift = get(hObj, 'Value');
        PlnShift = round(PlnShift);
        str = sprintf('Slice shift: %1d', PlnShift);
        set(hPlnShiftSliderText, 'String', str);
        
        UpdateDisplay;   
    end

    function continuous_pln_shift_slider_callback(hObj, eventdata)
        PlnShift = get(hObj, 'Value');
        PlnShift = round(PlnShift);
        str = sprintf('Slice shift: %1d', PlnShift);
        set(hPlnShiftSliderText, 'String', str);
        
        UpdateDisplay;   
    end

    function col_shift_slider_callback(hObj, eventdata)
        ColShift = get(hObj, 'Value');
        ColShift = round(ColShift);
        str = sprintf('Column shift: %1d', ColShift);
        set(hColShiftSliderText, 'String', str);
        
        UpdateDisplay;   
    end

    function continuous_col_shift_slider_callback(hObj, eventdata)
        ColShift = get(hObj, 'Value');
        ColShift = round(ColShift);
        str = sprintf('Column shift: %1d', ColShift);
        set(hColShiftSliderText, 'String', str);
        
        UpdateDisplay; 
    end

    function row_shift_slider_callback(hObj, eventdata)
        RowShift = get(hObj, 'Value');
        RowShift = round(RowShift);
        str = sprintf('Row shift: %1d', RowShift);
        set(hRowShiftSliderText, 'String', str);
        
        UpdateDisplay; 
    end

    function continuous_row_shift_slider_callback(hObj, eventdata)
        RowShift = get(hObj, 'Value');
        RowShift = round(RowShift);
        str = sprintf('Row shift: %1d', RowShift);
        set(hRowShiftSliderText, 'String', str);
        
        UpdateDisplay; 
    end

    % Deleting the figure also deletes its child objects and causes the function to return
    function pushbutton_callback(hObj, eventdata)
        CreateRegisteredArray;
        
        pause(0.25);
        
        delete(hf);
    end

    function my_closefcn(hObj, eventdata) 
        CreateRegisteredArray;
        
        pause(0.25);
        
        delete(hf);
    end

    function my_keypressfcn(hObj, eventdata)
        currChar = get(hf, 'CurrentKey');
        
        if (isequal(currChar, 'return') || isequal(currChar, 'escape'))   
            CreateRegisteredArray;
            
            pause(0.25);
            
            delete(hf);
        end
    end

    function CreateRegisteredArray
        Registered = padarray(Moving, [Pads, Pads, Pads], 0, 'both');
        Registered = circshift(Registered, [RowShift, ColShift, PlnShift]);
        Registered = Registered(1+Pads:end-Pads, 1+Pads:end-Pads, 1+Pads:end-Pads);
    end

    function UpdateDisplay
        F = Fixed(:, :, Pln);
        M = Moving(:, :, Pln - PlnShift);
        M = padarray(M, [Pads, Pads], 0, 'both');
        M = circshift(M, [RowShift, ColShift]);
        M = M(1+Pads:end-Pads, 1+Pads:end-Pads);
        imshowpair(F, M, 'falsecolor', 'Parent', hAxes12);
        
        F = Fixed(:, Col, :);
        F = squeeze(F);
        F = reshape(F, [ROWS, PLNS]);
        M = Moving(:, Col - ColShift, :);
        M = squeeze(M);
        M = reshape(M, [ROWS, PLNS]);
        M = padarray(M, [Pads, Pads], 0, 'both');
        M = circshift(M, [RowShift, PlnShift]);
        M = M(1+Pads:end-Pads, 1+Pads:end-Pads);
        imshowpair(F, M, 'falsecolor', 'Parent', hAxes13);   
        
        F = Fixed(Row, :, :);
        F = squeeze(F);
        F = reshape(F, [COLS, PLNS]);
        F = F';
        M = Moving(Row - RowShift, :, :);
        M = squeeze(M);
        M = reshape(M, [COLS, PLNS]);
        M = M';
        M = padarray(M, [Pads, Pads], 0, 'both');
        M = circshift(M, [PlnShift, ColShift]);
        M = M(1+Pads:end-Pads, 1+Pads:end-Pads);
        imshowpair(F, M, 'falsecolor', 'Parent', hAxes23);         
    end

end

