function [Imask] = extractMask(I, maskcolor)
% Create mask from red pixels

[ni, nj]=size(I(:,:,1));

Imask = zeros(ni, nj);

Imask(I(:,:,1) == maskcolor(1) & I(:,:,2) == maskcolor(2) & I(:,:,3) == maskcolor(3)) = 1;

end
