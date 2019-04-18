function Interpolation = pft_GetInterpolationType

% Initialize the output for possible immediate (default) return
Interpolation = 'linear';  

% Dimension the dialog box for the screen
x0 = 0.375;
y0 = 0.375;
wd = 0.25;
ht = 0.25;

% Create a main figure
hf = figure('Name', 'Processing decision', 'MenuBar', 'none', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [x0 y0 wd ht]);

% Create callbacks for exit by various means
set(hf, 'CloseRequestFcn', @my_closefcn);  % Trap the close button
set(hf, 'KeyPressFcn', @my_keypressfcn);   % Trap a RETURN or an ESCAPE from the keyboard

% Create a button group
hButtonGroup = uibuttongroup('Visible', 'off', 'Units', 'normalized', 'Position', [0.05 0.4 0.9 0.5], ...
                             'Title', 'Interpolation', 'FontName', 'Arial', 'FontUnits', 'points', ...
                             'FontSize', 16, 'Parent', hf);
        
uLinear  = uicontrol('Style', 'Radio', 'String', 'Linear (default)', 'Units', 'normalized', ...
                     'FontName', 'Arial', 'FontUnits', 'points', 'FontSize', 12, ...
                     'Position', [0.1 0.75 0.8 0.15], 'Parent', hButtonGroup, 'HandleVisibility', 'off');   
                  
uNearest = uicontrol('Style', 'Radio', 'String', 'Nearest-Neighbour', 'Units', 'normalized', ...
                     'FontName', 'Arial', 'FontUnits', 'points', 'FontSize', 12, ...
                     'Position', [0.1 0.50 0.8 0.15], 'Parent', hButtonGroup, 'HandleVisibility', 'off'); 
                 
uCubic   = uicontrol('Style', 'Radio', 'String', 'Cubic', 'Units', 'normalized', ...
                     'FontName', 'Arial', 'FontUnits', 'points', 'FontSize', 12, ...
                     'Position', [0.1 0.25 0.8 0.15], 'Parent', hButtonGroup, 'HandleVisibility', 'off');           
                 
set(hButtonGroup, 'SelectionChangeFcn', @InterpolationSelectionCallBack);

% Select and display an initial value
set(hButtonGroup, 'SelectedObject', uLinear); 

% Make the button group visible
set(hButtonGroup, 'Visible', 'on'); 
                                    
% Create an "APPLY" button
hApplyButton = uicontrol('Parent', hf, 'Style', 'Pushbutton', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.25], ...
                         'String', 'APPLY', 'FontName', 'Arial', 'FontUnits', 'points', 'FontSize', 18, 'FontWeight', 'bold', ...
                         'Callback', @ApplyButtonCallback);                                    
                               
% Block user interaction with MATLAB until the figure is destroyed and the function returns
uiwait(hf);

% Nested callbacks follow
function InterpolationSelectionCallBack(hObject, EventData)
  Code = get(EventData.NewValue, 'String');
  
  switch Code
    case 'Linear'
      Interpolation = 'linear';
    case 'Nearest-Neighbour'
      Interpolation = 'nearest';
    case 'Cubic'
      Interpolation = 'cubic';
  end
  
end

function ApplyButtonCallback(hObj, eventdata)
  delete(hf);
end

function my_closefcn(hObj, eventdata)
  delete(hf);
end

function my_keypressfcn(hObj, eventdata)
  currChar = get(hf, 'CurrentKey');
        
  if (isequal(currChar, 'return') || isequal(currChar, 'escape'))
    delete(hf);
  end
end

end







