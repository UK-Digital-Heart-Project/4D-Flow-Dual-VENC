function FileName = pft_NumberedFileName(n)

S = sprintf('%1d', n);
N = length(S);
P = 4 - N;
Z = repmat('0', [1, P]);

FileName = strcat('OFFLINE', '_', 'PRISMA', '_', Z, S, '.dcm');

end

