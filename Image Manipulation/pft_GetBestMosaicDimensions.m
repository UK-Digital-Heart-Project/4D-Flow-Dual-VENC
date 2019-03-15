function [ Rows, Cols ] = pft_GetBestMosaicDimensions(Wd, Ht, NP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to find suitable montage dimensions for an image stack with NP slices, each with width Wd and height Ht.                   %
%                                                                                                                                       %
% PFT - 18. 05. 2018.                                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the available space
ScreenSize = get(0, 'ScreenSize');

SW = ScreenSize(3); % Width
SH = ScreenSize(4); % Height

% Place limits on the number of tiles which can be created
Ht = double(Ht);
Wd = double(Wd);

MaxiRows = floor(SH/Ht);
MaxiCols = floor(SW/Wd);

% Construct a table of possible montage dimensions
[ CC, RR ] = meshgrid(1:MaxiCols, 1:MaxiRows);

Tiles = CC.*RR;

% Create tables to measure the 'footprint' and 'squareness' of the candidate montage dimensions
Summ = RR + CC;
Diff = abs(RR - CC);

% Initialise a table of possibilities
Possible = true([MaxiRows, MaxiCols]);

% Exclude the irrelevant entries from the 'tie-break' arrays
Exclude = (Tiles < NP);

Possible(Exclude) = false;
Summ(Exclude) = NaN;
Diff(Exclude) = NaN;

% Find the 'smallest' arrangement
MiniSumm = nanmin(Summ(:));

% Exclude possibilities which are larger
Exclude = (Summ > MiniSumm);

Possible(Exclude) = false;
Diff(Exclude) = NaN;

% Now locate the most 'square' option in the difference array
MiniDiff = nanmin(Diff(:));

% Exclude possibilities which are larger, as before
Possible(Diff > MiniDiff) = false;

% Search for the best remaining possibility by rows or columns, as appropriate, in order to keep the final montage compact
if (MaxiRows > MaxiCols)
  [ Rows, Cols ] = find(Possible == true, 1, 'first');
else
  Possible = Possible';  
    
  [ Cols, Rows ] = find(Possible == true, 1, 'first');
end

% Signal an error if no result is found
if isempty(Rows)
  h = warndlg('No suitable montage found.', 'Warning', 'modal');
  uiwait(h);
  delete(h);
end

end

