function [ Intercept, Slope ] = pft_GetVelocityScaling(Source)
% The Source may be a Dicom file or a Dicom header.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The resulting Intercept and Slope are double values.                                  %
% The usage is:                                                                         %
%               Velocity = Intercept + Slope*double(Grayscale)                          %
%                                                                                       %
% Pawel Tokarczuk - 25. 03. 2019.                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This snippet is used by Chris Rodgers in related code
if ischar(Source)
  Dictionary = dicomdict('get');
  
  Source = dicominfo(Source, 'Dictionary', Dictionary);
end 

% This is correct for processed images and accommodates an unsymmetrical display range
Intercept = double(Source.RescaleIntercept);
Slope     = double(Source.RescaleSlope);

% This applies to Siemens source images, whose Intercept and Slope are meaningless
TOL = 1.0e-2;

if (abs(Intercept + 4096.0) < TOL) || (abs(Slope - 2.0) < TOL)
  Venc = pft_GetVencFromHeader(Source);
  
  BS = Source.BitsStored;
  
  Intercept = - double(Venc);
  Slope     = 2.0*double(Venc)/double(2^BS);
end

end

