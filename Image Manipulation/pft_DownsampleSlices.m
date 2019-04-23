function Downsampled = pft_DownsampleSlices(Fine, SZ, TZ, Interpolation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover the sampled 3D image stack at the original slice locations after padding, up-sampling and co-registration.                        %
%                                                                                                                                           %
% PFT - 20. 04. 2018, modified - 23. 04. 2019.                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There is some padding that needs to be removed
E = 8;

% Determine the size and type of the i/p and o/p arrays
[ NR, NC, TP ] = size(Fine);

NP = numel(SZ);

Class = class(Fine);

% Create some FP working arrays
F = single(Fine);
D = zeros([NR, NC, NP], 'single'); 

% Now up-sample the i/p array one pixel at a time
parfor c = 1:NC
  for r = 1:NR
    V = F(r, c, :);
    V = reshape(V, 1, 1, TP);
    V = squeeze(V);
    W = interp1(TZ, V, SZ, Interpolation);
    W = reshape(W, 1, 1, NP);
    D(r, c, :) = W;
  end
end

% Now convert the FP downsampled array to the required class
Downsampled = eval(sprintf('%s(D)', Class));

Downsampled(:, :, 1:E)         = [];
Downsampled(:, :, end-E+1:end) = [];

end



