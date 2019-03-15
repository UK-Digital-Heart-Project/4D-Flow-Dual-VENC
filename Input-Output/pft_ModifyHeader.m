function Head = pft_ModifyHeader(Info, NewVenc, NewSeriesDescription, NewImageComments)

% First, fetch the private Siemens information
Head = Info;
Head = SiemensCsaParse(Head);

% Next, edit the sequence name to reflect the new Venc value
Head.SequenceName = strcat('FL', num2str(round(NewVenc)));

% Now add the new Venc to the private 'csa' field of the Siemens header
Head.csa.FlowVenc = NewVenc;

% And finally, some boilerplate options for co-registered or fused velocity images, including the all-important scaling slope and intercept
Head.BitsAllocated     = 16;
Head.BitsStored        = 16;
Head.HighBit           = 15;
Head.RescaleIntercept  = - NewVenc;
Head.RescaleSlope      = 2.0*NewVenc/double(2^16);
Head.RescaleType       = 'cm/s';
Head.SeriesDescription = NewSeriesDescription;
Head.ImageComments     = NewImageComments;

end







