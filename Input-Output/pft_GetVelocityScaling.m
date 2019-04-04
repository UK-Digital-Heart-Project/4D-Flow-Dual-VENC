function [ Intercept, Slope ] = pft_GetVelocityScaling(Source)
% The Source may be a Dicom file or a Dicom header.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The resulting Intercept and Slope are double values.                                  %
% The usage is:                                                                         %
%               Velocity = Intercept + Slope*double(Grayscale)                          %
%                                                                                       %
% Pawel Tokarczuk - 02. 04. 2019.                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This snippet is used by Chris Rodgers in related code
if ischar(Source)
  Dictionary = dicomdict('get');
  
  Source = dicominfo(Source, 'Dictionary', Dictionary);
end 

% The Series Description denotes a dual-Venc image if it contains the string 'RSS'
SD = Source.SeriesDescription;

if contains(SD, 'RSS')
  % This is correct for processed images and accommodates an unsymmetrical display range
  Intercept = double(Source.RescaleIntercept);
  Slope     = double(Source.RescaleSlope);
else
  % This applies to Siemens source images, whose Intercept and Slope are meaningless
  Venc = pft_GetVencFromHeader(Source);
  
  % This field is subject to transcription errors, and may be reset to 16, so disregard the Dicom value
  BS = 12;
  
  Intercept = - double(Venc);
  Slope     = 2.0*double(Venc)/double(2^BS);
end

end

