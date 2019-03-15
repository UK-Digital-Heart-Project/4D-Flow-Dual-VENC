function [ Data, Info ] = pft_ReadDicomCineStack(Folder)

% List the DICOM files in the nominated folder
Listing = dir(Folder);
Entries = { Listing.name };
Folders = [ Listing.isdir ];
Entries = Entries(~Folders);
Entries = sort(Entries);
Entries = Entries';

% Exit with an empty array if there are no DICOM files to find
if isempty(Entries)
  h = errordlg('No DICOM files found', 'Input error', 'modal');
  uiwait(h);
  delete(h);
  Data = [];
  Info = {};
  return;
end

% Proceed to normal execution
NFILES = int32(numel(Entries));

Info = cell(NFILES, 1);

% Sort the files by slice location and filename - this will need to be adapted later, when the inversion times TI become available
SliceLocations = zeros([NFILES, 1], 'double');
TriggerTimes   = zeros([NFILES, 1], 'double');

TOL = 1.0e-2;

wb = waitbar(0, 'Reading image headers ... ');

Dictionary = dicomdict('get');

for n = 1:NFILES
  Info{n} = dicominfo(fullfile(Folder, Entries{n}), 'Dictionary', Dictionary, 'UseVRHeuristic', false);
  
  SliceLocations(n) = TOL*round(Info{n}.SliceLocation/TOL);
  TriggerTimes(n)   = TOL*round(Info{n}.TriggerTime/TOL);
  
  waitbar(double(n)/double(NFILES), wb, sprintf('Read %1d of %1d headers.', n, NFILES));
end

waitbar(1, wb, 'Read all headers.');

delete(wb);

% Perform the sorting here
SortingTable = horzcat(SliceLocations, TriggerTimes);

[ SortedTable, Index ] = sortrows(SortingTable, [1, 2]);

Entries = Entries(Index);

Info = Info(Index);

% Determine how many times and slices there are
NSLICES = int32(numel(unique(SliceLocations)));
NEPOCHS = int32(idivide(NFILES, NSLICES));

% Exit in the event of an incomplete data set
if (mod(NFILES, NSLICES) ~= 0)
  h = errordlg('Incomplete data set', 'Input error', 'modal');
  uiwait(h);
  delete(h);
  Data = [];
  Info = {};
  return;
end

% Now read in the images
Sample = squeeze(dicomread(fullfile(Folder, Entries{1}), 'UseVRHeuristic', false));

[ NROWS, NCOLS ] = size(Sample);

Kind = class(Sample);

Data = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], Kind);

Count = 1;

wb = waitbar(0, 'Reading image data');

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Image = squeeze(dicomread(fullfile(Folder, Entries{Count}), 'UseVRHeuristic', false));
    Data(:, :, e, s) = Image;
    waitbar(double(Count)/double(NFILES), wb, sprintf('Read %1d of %1d images.', Count, NFILES));
    Count = Count + 1;
  end
end

waitbar(1, wb, 'Read all images.');

delete(wb);

% Finally, squeeze out any singleton dimensions - the image might be purely structural rather than dynamic
Data = squeeze(Data);

end
    
