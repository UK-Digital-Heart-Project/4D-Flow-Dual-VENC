function Venc = pft_GetVencFromHeader(Info)

% Fetch the private Siemens information and extract the value from there, rather than use the encoded sequence name
Head = Info;
Head = SiemensCsaParse(Head);

Venc = Head.csa.FlowVenc;

end

