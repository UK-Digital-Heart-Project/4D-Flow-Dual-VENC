function Mosaic = pft_MosaicImages(Stack, Rows, Cols, Wd, Ht)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a stack of images to a mosaic display.                                                                                      %
% The mosaic has dimensions of (Rows, Cols), with each tile having diemsnions of Wd and Ht, in an obvious notation.                 %
%                                                                                                                                   %
% PFT - 18. 05. 2018.                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fetch the dimensions of a 2-D slice or mosaic - the dimensions will always be multiples of 128
[ NR, NC, NP ] = size(Stack);   

% The case of a 1x1 image stack is trivial
if (NP == 1)
  Mosaic = Stack;  
  return;
end

% Initialise an o/p array
Class  = class(Stack);
Mosaic = zeros(Ht*Rows, Wd*Cols, Class);

% Write the planes of the stack to the appropriate tiles of the mosaic, making allowance for possible empty tiles at the end
P = 1;

InnerBreak = false;

for R = 1:Rows    
  for C = 1:Cols      
    ur = Ht*R;
    lr = ur - Ht + 1;
    uc = Wd*C;
    lc = uc - Wd + 1;
  
    Mosaic(lr:ur, lc:uc) = Stack(:, :, P);
    
    P = P + 1;
    
    if (P > NP)
      InnerBreak = true;
      break;
    end
  end
  
  if (InnerBreak == true)
    break;
  end   
end

end


  


