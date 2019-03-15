function [ Target, SZ, TZ ] = pft_InterpolateSlices(Source, DR, DZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to create a 3D array with isotropic pixels from a 2D stack, isotropic in-plane but thick in the slice direction.               %
%                                                                                                                                           %
% PFT - 20. 04. 2018.                                                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the i/p array size and type
[ NR, NC, NP ] = size(Source);

Class = class(Source);

% Create some padding in the slice direction to allow for out-of-slice rotations during co-registration
E = 8;

Roof = zeros([NR, NC, E], Class);
Deck = zeros([NR, NC, E], Class);

Source = cat(3, Deck, Source, Roof);

% Create an array of i/p z-coordinates, with zero at the first (original and un-padded) plane
SZ = DZ*(-E:NP-1+E);

% Now, an array of o/p z-coordinates, spaced at the in-plane resolution and long enough to cover the i/p image stack
LL = ceil(-E*double(DZ/DR));
UL = floor((NP-1+E)*double(DZ/DR));

TZ = DR*(LL:UL);
TP = numel(TZ);

% Create some FP working arrays
S = single(Source);
T = zeros([NR, NC, TP], 'single'); 

% Now up-sample the i/p array one pixel at a time
parfor c = 1:NC
  for r = 1:NR
    V = S(r, c, :);
    V = reshape(V, 1, 1, NP+2*E);
    V = squeeze(V);
    W = interp1(SZ, V, TZ, 'spline');
    W = reshape(W, 1, 1, TP);
    T(r, c, :) = W;
  end  
end

% Now convert the FP target array to the required class
Target = eval(sprintf('%s(T)', Class));

end



