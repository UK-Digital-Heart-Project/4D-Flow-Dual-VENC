%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the workspace

clear all
close all
clc

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate some input folders

if ispc
  Username = getenv('Username');
  Home = fullfile('C:', 'Users', Username, 'Desktop');
elseif isunix || ismac
  [ Status, CmdOut ] = system('whoami');
  Home = fullfile('home', CmdOut, 'Desktop');
end  
  
StartPath = uigetdir(Home, 'Select a top-level folder with scan folders inside');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominate a single output root folder

% 01. Sub-folders will be created in the worker function
MergedRoot = uigetdir(StartPath, 'Root folder for OUTPUT files');

if ~ischar(MergedRoot)
  h = msgbox('No folder chosen', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

%% Fetch the interpolation type

Interpolation = pft_GetInterpolationType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the worker function - this has been written so that it can be called programmatically for multiple acquisitions

pft_NonRigidThreeDimensionalMergingFunction(LoVencMagnitudeSource, HiVencMagnitudeSource, LoVencPhaseSource, HiVencPhaseSource, MergedRoot, Interpolation);