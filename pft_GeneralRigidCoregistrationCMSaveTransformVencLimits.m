%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the workspace
clear all
close all
clc

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate the source folders
StartPath = 'S:\4-D Flow MRI\Studies';

% 01. Low-Venc Magnitude
LoVencMagnitudeSource = uigetdir(StartPath, 'Low-Venc MAGNITUDE folder');

if ~ischar(LoVencMagnitudeSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end 

% 02. High-Venc Magnitude
HiVencMagnitudeSource = uigetdir(StartPath, 'High-Venc MAGNITUDE folder');

if ~ischar(HiVencMagnitudeSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end
  
% 03. Low-Venc Phase
LoVencPhaseSource = uigetdir(StartPath, 'Low-Venc PHASE folder');

if ~ischar(LoVencPhaseSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end 

% 04. High-Venc Phase
HiVencPhaseSource = uigetdir(StartPath, 'High-Venc PHASE folder');

if ~ischar(HiVencPhaseSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decide on the image interpolation:
%  Linear is the default;
%  Nearest-neighbour should avoid the wrapping-boundary artefact;
%  Cubic may give a smoother result than either of the other options, albeit with some residual wrapping-boundary artefact
Interpolation = pft_GetInterpolationType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate some o/p folders
Root = uigetdir(StartPath, 'Root folder for OUTPUT files');

if ~ischar(Root)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

% 01. Motion-corrected Low-Venc modulus
MoCoLoVencMagnitudeTarget = fullfile(Root, 'OFFLINE - MOCO LOW-VENC MAGNITUDE - RIGID');

if (exist(MoCoLoVencMagnitudeTarget, 'dir') ~= 7)
  mkdir(MoCoLoVencMagnitudeTarget);
end

% 02. Motion-corrected Low-Venc velocity
MoCoLoVencVelocityTarget = fullfile(Root, 'OFFLINE - MOCO LOW-VENC VELOCITY - RIGID');

if (exist(MoCoLoVencVelocityTarget, 'dir') ~= 7)
  mkdir(MoCoLoVencVelocityTarget);
end

% 03. Motion-corrected fused velocity
MoCoFusedVelocityTarget = fullfile(Root, 'OFFLINE - MOCO FUSED VELOCITY - RIGID');

if (exist(MoCoFusedVelocityTarget, 'dir') ~= 7)
  mkdir(MoCoFusedVelocityTarget);
end

% 04. Filtered motion-corrected fused velocity
FilteredMoCoFusedVelocityTarget = fullfile(Root, 'OFFLINE - FILTERED MOCO FUSED VELOCITY - RIGID');

if (exist(FilteredMoCoFusedVelocityTarget, 'dir') ~= 7)
  mkdir(FilteredMoCoFusedVelocityTarget);
end

% 05. Discrepancy between the motion-corrected fused velocity and the High-Venc velocity
DiscrepancyTarget = fullfile(Root, 'OFFLINE - MOCO FUSED VELOCITY DISCREPANCY - RIGID');

if (exist(DiscrepancyTarget, 'dir') ~= 7)
  mkdir(DiscrepancyTarget);
end

% 06. Twice-corrected fused velocity
TwiceCorrectedVelocityTarget = fullfile(Root, 'OFFLINE - TWICE-CORRECTED FUSED VELOCITY - RIGID');

if (exist(TwiceCorrectedVelocityTarget, 'dir') ~= 7)
  mkdir(TwiceCorrectedVelocityTarget);
end

% 07. Filtered twice-corrected fused velocity
FilteredTwiceCorrectedVelocityTarget = fullfile(Root, 'OFFLINE - FILTERED TWICE-CORRECTED FUSED VELOCITY - RIGID');

if (exist(FilteredTwiceCorrectedVelocityTarget, 'dir') ~= 7)
  mkdir(FilteredTwiceCorrectedVelocityTarget);
end

% 08. Residual between the twice-corrected fused velocity and the High-Venc velocity
ResidualTarget = fullfile(Root, 'OFFLINE - TWICE-CORRECTED FUSED VELOCITY RESIDUAL - RIGID');

if (exist(ResidualTarget, 'dir') ~= 7)
  mkdir(ResidualTarget);
end

% 09. Co-registration mosaic screenshots, both modulus and phase
ScreenshotTarget = fullfile(Root, 'OFFLINE - SCREENSHOTS - RIGID');

if (exist(ScreenshotTarget, 'dir') ~= 7)
  mkdir(ScreenshotTarget);
end

% 10. Manual shift corrections and transformation matrices
TransformationTarget = fullfile(Root, 'OFFLINE - TRANSFORMATION MATRICES - RIGID');

if (exist(TransformationTarget, 'dir') ~= 7)
  mkdir(TransformationTarget);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in the source data
[ LoVencMagnitudeData, LoVencMagnitudeInfo ] = pft_ReadDicomCineStack(LoVencMagnitudeSource);
[ HiVencMagnitudeData, HiVencMagnitudeInfo ] = pft_ReadDicomCineStack(HiVencMagnitudeSource);

[ LoVencPhaseData, LoVencPhaseInfo ] = pft_ReadDicomCineStack(LoVencPhaseSource);
[ HiVencPhaseData, HiVencPhaseInfo ] = pft_ReadDicomCineStack(HiVencPhaseSource);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract some basic dimensions from the High-Venc Magnitude image (and the first image slice in the stack)
[ NROWS, NCOLS, NEPOCHS, NSLICES ] = size(HiVencMagnitudeData);

Wd = NCOLS;
Ht = NROWS;
NP = NSLICES;

[ Rows, Cols ] = pft_GetBestMosaicDimensions(Wd, Ht, NP);

PS = HiVencMagnitudeInfo{1}.PixelSpacing(1);
PS = 0.01*round(PS/0.01);

DR = PS(1);

ST = HiVencMagnitudeInfo{1}.SliceThickness;
ST = 0.01*round(ST/0.01);

DZ = ST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert the phase images to velocity units (cm/s)
[ Intercept, Slope ] = pft_GetVelocityScaling(LoVencPhaseInfo{1});

LoVencVelocity = Intercept + Slope*double(LoVencPhaseData);

LoVenc = - Intercept;

[ Intercept, Slope ] = pft_GetVelocityScaling(HiVencPhaseInfo{1});

HiVencVelocity = Intercept + Slope*double(HiVencPhaseData);

HiVenc = - Intercept;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assign some o/p arrays
MoCoLoVencMagnitude = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');
MoCoLoVencVelocity  = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For each epoch, co-register the Low-Venc modulus image to the High-Venc modulus, then apply the same transformation to the Low-Venc phase image
iptsetpref('ImshowBorder', 'tight');

f = figure('Name', 'Co-registration', 'MenuBar', 'none', 'NumberTitle', 'off');
axis off, grid off, box off, hold on

[ Optimizer, Metric ] = imregconfig('Monomodal');

RowShift = 0;
ColShift = 0;
PlnShift = 0;

for e = 1:NEPOCHS
  % Interpolate the fixed (High-Venc) image to cubical voxels  
  Fixed = squeeze(HiVencMagnitudeData(:, :, e, :));  
  [ FixedFine, SZ, TZ ] = pft_InterpolateSlices(Fixed, DR, DZ, Interpolation);
  
  % Interpolate the moving (Low-Venc) image to cubical voxels
  Moving = squeeze(LoVencMagnitudeData(:, :, e, :));
  [ MovingFine, SZ, TZ ] = pft_InterpolateSlices(Moving, DR, DZ, Interpolation);
  
  % Save 2 modulus image mosaics for later display
  AA = pft_MosaicImages(Fixed, Rows, Cols, Wd, Ht);
  BB = pft_MosaicImages(Moving, Rows, Cols, Wd, Ht);
  
  % Co-register the modulus images approximately and manually - apply the correction to the moving image in-place
  [ MovingFine, RowShift, ColShift, PlnShift ] = pft_CoregisterInteractively(MovingFine, FixedFine, RowShift, ColShift, PlnShift);
  
  % Save the manual shift corrections in a CSV file
  Leaf = sprintf('Manual-Shift-%03d.csv', e);
  
  T = table(RowShift, ColShift, PlnShift, 'VariableNames', { 'Row_Shift', 'Col_Shift', 'Pln_Shift' });
  
  writetable(T, fullfile(TransformationTarget, Leaf));  
  
  % Now perform the rigid co-registration and save the transformation for re-use
  Transform = imregtform(MovingFine, FixedFine, 'rigid', Optimizer, Metric);
  
  [ RegisteredMovingFine, DisposableReferenceObject ] = imwarp(MovingFine, Transform, Interpolation, 'OutputView', imref3d(size(MovingFine)));
    
  % Downsample the result in the z-direction  
  RegisteredMoving = pft_DownsampleSlices(RegisteredMovingFine, SZ, TZ, Interpolation);
  
  % Save the transformation to a MAT file for auditing and possible later re-use (outside this script)
  Leaf = sprintf('Transform-%03d.mat', e);
  
  save(fullfile(TransformationTarget, Leaf), 'Transform');
  
  MoCoLoVencMagnitude(:, :, e, :) = RegisteredMoving;
  
  % Create another modulus image mosaic for later display
  CC = pft_MosaicImages(RegisteredMoving, Rows, Cols, Wd, Ht);
  
  % Display the original High-Venc and Low-Venc modulus image mosaics together
  imshowpair(AA, BB, 'falsecolor');
  pause(0.25);
  imshowpair(AA, BB, 'falsecolor');
  text(8, 8, 'Hi-Venc and Lo-Venc Modulus', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 24, sprintf('Epoch %1d', e), 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  
  % Save the pink-and-green overlay for later review
  FilePath = fullfile(ScreenshotTarget, sprintf('Hi-Venc and Lo-Venc Modulus Epoch %02d.png', e));
  F = getframe(f);
  X = F.cdata;
  imwrite(X, FilePath);
  
  pause(2.0);
  
  % Display the original High-Venc and the motion-corrected Low-Venc modulus image mosaics together
  imshowpair(AA, CC, 'falsecolor');
  pause(0.25);
  imshowpair(AA, CC, 'falsecolor');
  text(8, 8, 'Hi-Venc and Lo-Venc Modulus', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 24, sprintf('Epoch %1d', e), 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 40, 'Lo-Venc Motion-Corrected', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  
  % Save the pink-and-green overlay for later review
  FilePath = fullfile(ScreenshotTarget, sprintf('Hi-Venc and MoCo Lo-Venc Modulus Epoch %02d.png', e));
  F = getframe(f);
  X = F.cdata;
  imwrite(X, FilePath);
  
  pause(2.0);
  
  % Interpolate the Low-Venc velocity image to cubical voxels
  Moving = squeeze(LoVencVelocity(:, :, e, :));
  [ MovingFine, SZ, TZ ] = pft_InterpolateSlices(Moving, DR, DZ, Interpolation);
  
  % Apply the previously applied manual circular shift (with padding beforehand, and trimming afterwards)
  PADS = 32;
  
  MovingFine = padarray(MovingFine, [ PADS, PADS, PADS ], 0, 'both');
  MovingFine = circshift(MovingFine, [ RowShift, ColShift, PlnShift ]);
  MovingFine = MovingFine(1+PADS:end-PADS, 1+PADS:end-PADS, 1+PADS:end-PADS); 
  
  % Apply the previously computed rigid transformation
  [ RegisteredMovingFine, DisposableReferenceObject ] = imwarp(MovingFine, Transform, Interpolation, 'OutputView', imref3d(size(MovingFine)));
  
  % Downsample the result in the z-direction
  RegisteredMoving = pft_DownsampleSlices(RegisteredMovingFine, SZ, TZ, Interpolation);
  
  % Extract an epoch from the original High-Venc cine-stack
  Fixed = squeeze(HiVencVelocity(:, :, e, :));
  
  % Save 2 velocity image mosaics for later display
  AA = pft_MosaicImages(Fixed, Rows, Cols, Wd, Ht);
  BB = pft_MosaicImages(Moving, Rows, Cols, Wd, Ht);
  
  % Insert an epoch of the motion-corrected velocity image into the output cine-stack
  MoCoLoVencVelocity(:, :, e, :) = RegisteredMoving;
  
  % Create another velocity image mosaic for later display
  CC = pft_MosaicImages(RegisteredMoving, Rows, Cols, Wd, Ht);
  
  % Display the original High-Venc and Low-Venc velocity image mosaics together
  imshowpair(AA, BB, 'falsecolor');
  pause(0.25);
  imshowpair(AA, BB, 'falsecolor');
  text(8, 8, 'Hi-Venc and Lo-Venc Velocity', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 24, sprintf('Epoch %1d', e), 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  
  % Save the pink-and-green overlay for later review
  FilePath = fullfile(ScreenshotTarget, sprintf('Hi-Venc and Lo-Venc Velocity Epoch %02d.png', e));
  F = getframe(f);
  X = F.cdata;
  imwrite(X, FilePath);
  
  pause(2.0);
  
  % Display the original High-Venc and the motion-corrected Low-Venc velocity image mosaics together
  imshowpair(AA, CC, 'falsecolor');
  pause(0.25);
  imshowpair(AA, CC, 'falsecolor');
  text(8, 8, 'Hi-Venc and Lo-Venc Velocity', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 24, sprintf('Epoch %1d', e), 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  text(8, 40, 'Lo-Venc Motion-Corrected', 'Color', [1 1 0], 'FontSize', 12, 'FontWeight', 'bold');
  
  % Save the pink-and-green overlay for later review
  FilePath = fullfile(ScreenshotTarget, sprintf('Hi-Venc and MoCo Lo-Venc Velocity Epoch %02d.png', e));
  F = getframe(f);
  X = F.cdata;
  imwrite(X, FilePath);   
  
  pause(2.0);
end

delete(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create the fused velocity images - begin with an estimate of the wrapping correction, limited at either end to 4 times 2*LoVenc
Delta = 2.0*LoVenc*round(0.5*(HiVencVelocity - MoCoLoVencVelocity)/LoVenc);

Exclude = (abs(Delta) > 8.0*LoVenc);    
Include = ~Exclude;

% Rescale the o/p to 16-bits
MoCoFusedVelocity = MoCoLoVencVelocity;

MoCoFusedVelocity(Include) = MoCoLoVencVelocity(Include) + Delta(Include);

LL = LoVenc*floor(min(MoCoFusedVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(MoCoFusedVelocity(:))/LoVenc);

NewFusedVenc = max(abs(LL), abs(UL));

FusedPhaseData = uint16(double(2^15)*(MoCoFusedVelocity/NewFusedVenc + 1.0));

% Create the discrepancy image
DiscrepantVelocity = MoCoFusedVelocity - HiVencVelocity;

LL = LoVenc*floor(min(DiscrepantVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(DiscrepantVelocity(:))/LoVenc);

DiscrepantNewVenc = max(abs(LL), abs(UL));

DiscrepantPhaseData = uint16(double(2^15)*(DiscrepantVelocity/DiscrepantNewVenc + 1.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create the FILTERED fused velocity images
FilteredMoCoFusedVelocity = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');

Kernel = [ 1.0 2.0 1.0 ; ...
           2.0 4.0 2.0 ; ...
           1.0 2.0 1.0 ]/16.0;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Part = MoCoFusedVelocity(:, :, e, s);
    
    FilteredMoCoFusedVelocity(:, :, e, s) = conv2(Part, Kernel, 'same');
  end
end

LL = LoVenc*floor(min(FilteredMoCoFusedVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(FilteredMoCoFusedVelocity(:))/LoVenc);

NewFilteredFusedVenc = max(abs(LL), abs(UL));

FilteredFusedPhaseData = uint16(double(2^15)*(FilteredMoCoFusedVelocity/NewFilteredFusedVenc + 1.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Create the TWICE-CORRECTED fused velocity images to account for interpolation effects at wrapping boundaries after co-registration
Epsilon = LoVenc*round((HiVencVelocity - MoCoFusedVelocity)/LoVenc);

Exclude = (abs(Epsilon) > LoVenc);    
Include = ~Exclude;

% Rescale the o/p to 16-bits
TwiceCorrectedFusedVelocity = MoCoFusedVelocity;

TwiceCorrectedFusedVelocity(Include) = MoCoFusedVelocity(Include) + Epsilon(Include);

LL = LoVenc*floor(min(TwiceCorrectedFusedVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(TwiceCorrectedFusedVelocity(:))/LoVenc);

TwiceCorrectedNewFusedVenc = max(abs(LL), abs(UL));

TwiceCorrectedFusedPhaseData = uint16(double(2^15)*(TwiceCorrectedFusedVelocity/TwiceCorrectedNewFusedVenc + 1.0));

% Create the residual image
ResidualVelocity = TwiceCorrectedFusedVelocity - HiVencVelocity;

LL = LoVenc*floor(min(ResidualVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(ResidualVelocity(:))/LoVenc);

ResidualNewVenc = max(abs(LL), abs(UL));

ResidualPhaseData = uint16(double(2^15)*(ResidualVelocity/ResidualNewVenc + 1.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create the FILTERED twice-corrected velocity images
FilteredTwiceCorrectedFusedVelocity = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');

Kernel = [ 1.0 2.0 1.0 ; ...
           2.0 4.0 2.0 ; ...
           1.0 2.0 1.0 ]/16.0;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Part = TwiceCorrectedFusedVelocity(:, :, e, s);
    
    FilteredTwiceCorrectedFusedVelocity(:, :, e, s) = conv2(Part, Kernel, 'same');
  end
end

LL = LoVenc*floor(min(FilteredTwiceCorrectedFusedVelocity(:))/LoVenc);
UL = LoVenc*ceil(max(FilteredTwiceCorrectedFusedVelocity(:))/LoVenc);

FilteredTwiceCorrectedNewFusedVenc = max(abs(LL), abs(UL));

FilteredTwiceCorrectedFusedPhaseData = uint16(double(2^15)*(FilteredTwiceCorrectedFusedVelocity/FilteredTwiceCorrectedNewFusedVenc + 1.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 01. Write out the motion-corrected Lo-Venc modulus images
Dictionary = dicomdict('get');

wb = waitbar(0, 'Writing out MoCo Lo-Venc modulus images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = LoVencMagnitudeInfo{n};
  
    Head.RescaleType       = 'Grayscale';
    Head.SeriesDescription = 'Synthetic RSS image';
    Head.ImageComments     = 'MoCo Lo-Venc Magnitude image';
  
    OutputPathName = fullfile(MoCoLoVencMagnitudeTarget, pft_NumberedFileName(n));
  
    dicomwrite(uint16(MoCoLoVencMagnitude(:, :, e, s)), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 02. Write out the motion-corrected Lo-Venc velocity images
MoCoLoVencPhaseData = uint16(double(2^15)*(MoCoLoVencVelocity/LoVenc + 1.0));

wb = waitbar(0, 'Writing out MoCo Lo-Venc velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, LoVenc, 'Synthetic RSS image', '16-bit MoCo Low-Venc phase image');
    
    OutputPathName = fullfile(MoCoLoVencVelocityTarget, pft_NumberedFileName(n));
  
    dicomwrite(MoCoLoVencPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 03. Write out the motion-corrected fused velocity images
wb = waitbar(0, 'Writing out fused velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, NewFusedVenc, 'Synthetic RSS image', '16-bit fused phase image');
  
    OutputPathName = fullfile(MoCoFusedVelocityTarget, pft_NumberedFileName(n));
  
    dicomwrite(FusedPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 04. Write out the filtered motion-corrected fused velocity
wb = waitbar(0, 'Writing out filtered fused velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, NewFilteredFusedVenc, 'Synthetic RSS image', '16-bit filtered fused phase image');
  
    OutputPathName = fullfile(FilteredMoCoFusedVelocityTarget, pft_NumberedFileName(n));
  
    dicomwrite(FilteredFusedPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 05. Write out the discrepant velocity images
wb = waitbar(0, 'Writing out discrepant velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, DiscrepantNewVenc, 'Synthetic RSS image', '16-bit discrepant phase image');
    
    OutputPathName = fullfile(DiscrepancyTarget, pft_NumberedFileName(n));
  
    dicomwrite(DiscrepantPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 06. Write out the twice-corrected fused velocity images
wb = waitbar(0, 'Writing out twice-corrected fused velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, TwiceCorrectedNewFusedVenc, 'Synthetic RSS image', '16-bit twice-corrected fused phase image');
  
    OutputPathName = fullfile(TwiceCorrectedVelocityTarget, pft_NumberedFileName(n));
  
    dicomwrite(TwiceCorrectedFusedPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 07. Write out the filtered twice-corrected fused velocity images
wb = waitbar(0, 'Writing out filtered twice-corrected fused velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, FilteredTwiceCorrectedNewFusedVenc, 'Synthetic RSS image', '16-bit filtered twice-corrected fused phase image');
  
    OutputPathName = fullfile(FilteredTwiceCorrectedVelocityTarget, pft_NumberedFileName(n));
  
    dicomwrite(FilteredTwiceCorrectedFusedPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 08. Write out the residual velocity images
wb = waitbar(0, 'Writing out residual velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, ResidualNewVenc, 'Synthetic RSS image', '16-bit residual phase image');
    
    OutputPathName = fullfile(ResidualTarget, pft_NumberedFileName(n));
  
    dicomwrite(ResidualPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 09. Screenshots have already been saved on-the-fly during the co-registration phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 10. Manual shifts and transformation matrices have also been saved already

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write out a small summary file
fid = fopen(fullfile(Root, 'Summary - Rigid Co-Registration.txt'), 'wt');

fprintf(fid, 'Lo-Venc Magnitude source folder: %s\n', LoVencMagnitudeSource);
fprintf(fid, 'Hi-Venc Magnitude source folder: %s\n', HiVencMagnitudeSource);
fprintf(fid, 'Lo-Venc Phase source folder:     %s\n', LoVencPhaseSource);
fprintf(fid, 'Hi-Venc Phase source folder:     %s\n', HiVencPhaseSource);

fprintf(fid, '\n');

fprintf(fid, 'Output root folder: %s\n', Root);

fprintf(fid, '\n');

fprintf(fid, 'Original Low Venc  = %.2f cm/s\n', LoVenc);
fprintf(fid, 'Original High Venc = %.2f cm/s\n', HiVenc);

fprintf(fid, '\n');

fprintf(fid, 'Synthetic Once-Corrected Fused Venc  = %.2f cm/s\n', NewFusedVenc);
fprintf(fid, 'Filtered Once-Corrected Fused Venc   = %.2f cm/s\n', NewFilteredFusedVenc);
fprintf(fid, 'Synthetic Twice-Corrected Fused Venc = %.2f cm/s\n', TwiceCorrectedNewFusedVenc);
fprintf(fid, 'Filtered Twice-Corrected Fused Venc  = %.2f cm/s\n', FilteredTwiceCorrectedNewFusedVenc);

fprintf(fid, '\n');

fprintf(fid, 'Once-Corrected Discrepant Venc = %.2f cm/s\n', DiscrepantNewVenc);
fprintf(fid, 'Twice-Corrected Residual Venc  = %.2f cm/s\n', ResidualNewVenc);

fprintf(fid, '\n');

fprintf(fid, 'Co-registration:                           Rigid\n');
fprintf(fid, 'Interpolation used during co-registration: %s\n', Interpolation);

fprintf(fid, '\n');

fprintf(fid, 'For the 16-bit fused velocity image:\n');
fprintf(fid, '\n');
fprintf(fid, 'Offset = %1d\n', round(-NewFusedVenc));
fprintf(fid, 'Slope  = %.12f\n', 2.0*NewFusedVenc/double(2^16));

fprintf(fid, '\n');

fprintf(fid, 'For the 16-bit twice-corrected velocity image:\n');
fprintf(fid, '\n');
fprintf(fid, 'Offset = %1d\n', round(-TwiceCorrectedNewFusedVenc));
fprintf(fid, 'Slope  = %.12f\n', 2.0*TwiceCorrectedNewFusedVenc/double(2^16));

fclose(fid);

%% Write out a tidy XLSX file with information grouped into several tabs

FolderData = { 'INPUTS',             ' '; ...
               'Lo-Venc Magnitude',  LoVencMagnitudeSource; ...
               'Hi-Venc Magnitude',  HiVencMagnitudeSource; ...
               'Lo-Venc Phase',      LoVencPhaseSource; ...
               'Hi-Venc Phase',      HiVencPhaseSource; ...
               'OUTPUTS ',           ' ';
               'Output Root Folder', Root };
           
xlswrite(fullfile(Root, 'Processing Summary.xlsx'), FolderData, 'Data Folders');

ProcessingData = { 'Co-Registration', 'Interpolation'; ...
                   'Rigid',            Interpolation };
               
xlswrite(fullfile(Root, 'Processing Summary.xlsx'), ProcessingData, 'Processing');

VencData = { 'Image',                         'Venc [cm/s]',                      'Intercept',                          'Slope'; ...
             'Original Low-Venc',              LoVenc,                             - LoVenc,                             2.0*LoVenc/double(2^12); ...
             'Original High-Venc',             HiVenc,                             - HiVenc,                             2.0*HiVenc/double(2^12); ...
             'Once-Corrected Fused',           NewFusedVenc,                       - NewFusedVenc,                       2.0*NewFusedVenc/double(2^16); ...
             'Filtered Once-Corrected Fused',  NewFilteredFusedVenc,               - NewFilteredFusedVenc,               2.0*NewFilteredFusedVenc/double(2^16); ...
             'Twice-Corrected Fused',          TwiceCorrectedNewFusedVenc,         - TwiceCorrectedNewFusedVenc,         2.0*TwiceCorrectedNewFusedVenc/double(2^16); ...
             'Filtered Twice-Corrected Fused', FilteredTwiceCorrectedNewFusedVenc, - FilteredTwiceCorrectedNewFusedVenc, 2.0*FilteredTwiceCorrectedNewFusedVenc/double(2^16); ...
             'Once-Corrected Discrepancy',     DiscrepantNewVenc,                  - DiscrepantNewVenc,                  2.0*DiscrepantNewVenc/double(2^16); ...
             'Twice-Corrected Discrepancy',    ResidualNewVenc,                    - ResidualNewVenc,                    2.0*ResidualNewVenc/double(2^16) };
         
xlswrite(fullfile(Root, 'Processing Summary.xlsx'), VencData, 'Vencs and Velocity Scaling');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Signal completion
h = msgbox('All done !', 'Exit', 'modal');
uiwait(h);
delete(h);




