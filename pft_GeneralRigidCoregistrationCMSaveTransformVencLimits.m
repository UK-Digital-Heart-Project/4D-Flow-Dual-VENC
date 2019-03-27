%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the workspace
clear all
close all
clc

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate the source folders
StartPath = 'S:\4-D Flow MRI\Studies';

LoVencMagnSource = uigetdir(StartPath, 'Lo-Venc MAGNITUDE folder');

if ~ischar(LoVencMagnSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end 

HiVencMagnSource = uigetdir(StartPath, 'Hi-Venc MAGNITUDE folder');

if ~ischar(HiVencMagnSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end
  
LoVencPhaseSource = uigetdir(StartPath, 'Lo-Venc PHASE folder');

if ~ischar(LoVencPhaseSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end 

HiVencPhaseSource = uigetdir(StartPath, 'Hi-Venc PHASE folder');

if ~ischar(HiVencPhaseSource)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate some o/p folders
Root = uigetdir(StartPath, 'Root folder for OUTPUT files');

if ~ischar(Root)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

% Motion-corrected Low-Venc modulus
MoCoLoVencMagnTarget = fullfile(Root, 'OFFLINE - MOCO LOW-VENC MAGNITUDE - RIGID');

if (exist(MoCoLoVencMagnTarget, 'dir') ~= 7)
  mkdir(MoCoLoVencMagnTarget);
end

% Motion-corrected Low-Venc velocity
MoCoLoVencVelyTarget = fullfile(Root, 'OFFLINE - MOCO LOW-VENC VELOCITY - RIGID');

if (exist(MoCoLoVencVelyTarget, 'dir') ~= 7)
  mkdir(MoCoLoVencVelyTarget);
end

% Motion-corrected fused velocity
MoCoFusionVelyTarget = fullfile(Root, 'OFFLINE - MOCO FUSED VELOCITY - RIGID');

if (exist(MoCoFusionVelyTarget, 'dir') ~= 7)
  mkdir(MoCoFusionVelyTarget);
end

% Discrepancy between the motion-corrected fused velocity and the High-Venc velocity
DiscrepancyTarget = fullfile(Root, 'OFFLINE - MOCO FUSED VELOCITY DISCREPANCY - RIGID');

if (exist(DiscrepancyTarget, 'dir') ~= 7)
  mkdir(DiscrepancyTarget);
end

% Co-registration mosaic screenshots, both modulus and phase
ScreenshotTarget = fullfile(Root, 'OFFLINE - SCREENSHOTS - RIGID');

if (exist(ScreenshotTarget, 'dir') ~= 7)
  mkdir(ScreenshotTarget);
end

% Manual shift corrections and displacement fields
TransformTarget = fullfile(Root, 'TRANSFORMATION MATRICES');

if (exist(TransformTarget, 'dir') ~= 7)
  mkdir(TransformTarget);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in the source data
[ LoVencMagnData, LoVencMagnInfo ] = pft_ReadDicomCineStack(LoVencMagnSource);
[ HiVencMagnData, HiVencMagnInfo ] = pft_ReadDicomCineStack(HiVencMagnSource);

[ LoVencPhaseData, LoVencPhaseInfo ] = pft_ReadDicomCineStack(LoVencPhaseSource);
[ HiVencPhaseData, HiVencPhaseInfo ] = pft_ReadDicomCineStack(HiVencPhaseSource);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract some basic dimensions from the High-Venc magnitude image (and the first image slice in the stack)
[ NROWS, NCOLS, NEPOCHS, NSLICES ] = size(HiVencMagnData);

Wd = NCOLS;
Ht = NROWS;
NP = NSLICES;

[ Rows, Cols ] = pft_GetBestMosaicDimensions(Wd, Ht, NP);

Rows = Rows - 1;
Cols = Cols + 1;

PS = HiVencMagnInfo{1}.PixelSpacing(1);
PS = 0.01*round(PS/0.01);

DR = PS(1);

ST = HiVencMagnInfo{1}.SliceThickness;
ST = 0.01*round(ST/0.01);

DZ = ST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert the phase images to velocity units (cm/s)
HiVenc = pft_GetVencFromHeader(HiVencPhaseInfo{1});
LoVenc = pft_GetVencFromHeader(LoVencPhaseInfo{1});

HVBS = HiVencPhaseInfo{1}.BitsStored;
LVBS = LoVencPhaseInfo{1}.BitsStored;

HiVencVelocity = double(HiVenc)*(2.0*double(HiVencPhaseData)/double(2^HVBS) - 1.0);
LoVencVelocity = double(LoVenc)*(2.0*double(LoVencPhaseData)/double(2^LVBS) - 1.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assign some o/p arrays
MoCoLoVencMagn = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');
MoCoLoVencVely = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');

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
  Fixed = squeeze(HiVencMagnData(:, :, e, :));  
  [ FixedFine, SZ, TZ ] = pft_InterpolateSlices(Fixed, DR, DZ);
  
  % Interpolate the moving (Low-Venc) image to cubical voxels
  Moving = squeeze(LoVencMagnData(:, :, e, :));
  [ MovingFine, SZ, TZ ] = pft_InterpolateSlices(Moving, DR, DZ);
  
  % Save 2 modulus image mosaics for later display
  AA = pft_MosaicImages(Fixed, Rows, Cols, Wd, Ht);
  BB = pft_MosaicImages(Moving, Rows, Cols, Wd, Ht);
  
  % Co-register the modulus images approximately and manually - apply the correction to the moving image in-place
  [ MovingFine, RowShift, ColShift, PlnShift ] = pft_CoregisterInteractively(MovingFine, FixedFine, RowShift, ColShift, PlnShift);
  
  % Save the manual shift corrections in a CSV file
  Leaf = sprintf('Manual-Shift-%03d.csv', e);
  
  T = table(RowShift, ColShift, PlnShift, 'VariableNames', { 'Row_Shift', 'Col_Shift', 'Pln_Shift' });
  
  writetable(T, fullfile(TransformTarget, Leaf));  
  
  % Now perform the rigid co-registration and save the transformation for re-use
  Transform = imregtform(MovingFine, FixedFine, 'rigid', Optimizer, Metric);
  
  RegisteredMovingFine = imwarp(MovingFine, Transform, 'OutputView', imref3d(size(MovingFine)));
  
  % Downsample the result in the z-direction  
  RegisteredMoving = pft_DownsampleSlices(RegisteredMovingFine, SZ, TZ);
  
  % Save the transformation to a MAT file for auditing and possible later re-use (outside this script)
  Leaf = sprintf('Transform-%03d.mat', e);
  
  save(fullfile(TransformTarget, Leaf), 'Transform');
  
  MoCoLoVencMagn(:, :, e, :) = RegisteredMoving;
  
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
  [ MovingFine, SZ, TZ ] = pft_InterpolateSlices(Moving, DR, DZ);
  
  % Apply the previously applied manual circular shift (with padding beforehand, and trimming afterwards)
  PADS = 32;
  
  MovingFine = padarray(MovingFine, [ PADS, PADS, PADS ], 0, 'both');
  MovingFine = circshift(MovingFine, [ RowShift, ColShift, PlnShift ]);
  MovingFine = MovingFine(1+PADS:end-PADS, 1+PADS:end-PADS, 1+PADS:end-PADS); 
  
  % Apply the previously computed rigid transformation
  RegisteredMovingFine = imwarp(MovingFine, Transform, 'OutputView', imref3d(size(MovingFine)));
  
  % Downsample the result in the z-direction
  RegisteredMoving = pft_DownsampleSlices(RegisteredMovingFine, SZ, TZ);
  
  % Extract an epoch from the original High-Venc cine-stack
  Fixed = squeeze(HiVencVelocity(:, :, e, :));
  
  % Save 2 velocity image mosaics for later display
  AA = pft_MosaicImages(Fixed, Rows, Cols, Wd, Ht);
  BB = pft_MosaicImages(Moving, Rows, Cols, Wd, Ht);
  
  % Insert an epoch of the motion-corrected velocity image into the output cine-stack
  MoCoLoVencVely(:, :, e, :) = RegisteredMoving;
  
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
Delta = 2.0*LoVenc*round(0.5*(HiVencVelocity - MoCoLoVencVely)/LoVenc);

Exclude = (abs(Delta) > 8.0*LoVenc);    
Include = ~Exclude;

% Rescale the o/p to 16-bits
MoCoFusionVely = zeros([NROWS, NCOLS, NEPOCHS, NSLICES], 'double');

MoCoFusionVely(Include) = MoCoLoVencVely(Include) + Delta(Include);

LL = LoVenc*floor(min(MoCoFusionVely(:))/LoVenc);
UL = LoVenc*ceil(max(MoCoFusionVely(:))/LoVenc);

NewFusedVenc = max(abs(LL), abs(UL));

FusedPhaseData = uint16(double(2^15)*(MoCoFusionVely/NewFusedVenc + 1.0));

% Create the discrepancy image
DiscrepantVely = MoCoFusionVely - HiVencVelocity;

LL = LoVenc*floor(min(DiscrepantVely(:))/LoVenc);
UL = LoVenc*ceil(max(DiscrepantVely(:))/LoVenc);

DiscrepantVenc = max(abs(LL), abs(UL));

DiscrepantPhaseData = uint16(double(2^15)*(DiscrepantVely/DiscrepantVenc + 1.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write out the fused phase images first
Dictionary = dicomdict('get');

wb = waitbar(0, 'Writing out fused velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, NewFusedVenc, 'Synthetic RSS image', '16-bit fused phase image');
  
    OutputPathName = fullfile(MoCoFusionVelyTarget, pft_NumberedFileName(n));
  
    dicomwrite(FusedPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write out the motion-corrected Lo-Venc velocity images
MoCoLoVencPhaseData = uint16(double(2^15)*(MoCoLoVencVely/LoVenc + 1.0));

wb = waitbar(0, 'Writing out MoCo Lo-Venc velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, LoVenc, 'Synthetic RSS image', '16-bit MoCo Low-Venc phase image');
    
    OutputPathName = fullfile(MoCoLoVencVelyTarget, pft_NumberedFileName(n));
  
    dicomwrite(MoCoLoVencPhaseData(:, :, e, s), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write out the motion-corrected Lo-Venc modulus images
wb = waitbar(0, 'Writing out MoCo Lo-Venc modulus images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = LoVencMagnInfo{n};
  
    Head.RescaleType       = 'Grayscale';
    Head.SeriesDescription = 'Synthetic RSS image';
    Head.ImageComments     = 'MoCo Lo-Venc magnitude image';
  
    OutputPathName = fullfile(MoCoLoVencMagnTarget, pft_NumberedFileName(n));
  
    dicomwrite(uint16(MoCoLoVencMagn(:, :, e, s)), OutputPathName, Head, 'CreateMode', 'copy', 'Dictionary', Dictionary, 'WritePrivate', true);
  
    waitbar(double(n)/double(NFILES), wb, sprintf('%1d of %1d files written', n, NFILES));
    
    n = n + 1;
  end
end

waitbar(1, wb, sprintf('%1d of %1d files written', NFILES, NFILES));
pause(1.0);
delete(wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write out the discrepant velocity images
wb = waitbar(0, 'Writing out discrepant velocity images');

NFILES = NEPOCHS*NSLICES;

n = 1;

for s = 1:NSLICES
  for e = 1:NEPOCHS
    Head = pft_ModifyHeader(LoVencPhaseInfo{n}, DiscrepantVenc, 'Synthetic RSS image', '16-bit discrepant phase image');
    
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

%% Write out a small summary file
fid = fopen(fullfile(Root, 'Summary - Rigid Co-Registration.txt'), 'wt');

fprintf(fid, 'Lo-Venc Magnitude source folder: %s\n', LoVencMagnSource);
fprintf(fid, 'Hi-Venc Magnitude source folder: %s\n', HiVencMagnSource);
fprintf(fid, 'Lo-Venc Phase source folder:     %s\n', LoVencPhaseSource);
fprintf(fid, 'Hi-Venc Phase source folder:     %s\n', HiVencPhaseSource);

fprintf(fid, '\n');

fprintf(fid, 'Output root folder:              %s\n', Root);

fprintf(fid, '\n');

fprintf(fid, 'Low Venc   = %.2f cm/s\n', LoVenc);
fprintf(fid, 'High Venc  = %.2f cm/s\n', HiVenc);
fprintf(fid, 'Fused Venc = %.2f cm/s\n', NewFusedVenc);

fprintf(fid, '\n');

fprintf(fid, 'Residual Venc = %.2f cm/s\n', DiscrepantVenc);

fprintf(fid, '\n');

fprintf(fid, 'For the 16-bit fused velocity image:\n');
fprintf(fid, '\n');
fprintf(fid, 'Offset = %1d\n', round(-NewFusedVenc));
fprintf(fid, 'Slope  = %.12f\n', 2.0*NewFusedVenc/double(2^16));

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Signal completion
h = msgbox('All done !', 'Exit', 'modal');
uiwait(h);
delete(h);




