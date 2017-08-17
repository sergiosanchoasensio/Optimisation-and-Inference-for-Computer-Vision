% Extra Images
clearvars;

% name= 'ExtraImages/Lenna_20_toRestore.png';
% name= 'ExtraImages/Lenna_60_toRestore.png';
name= 'ExtraImages/Lenna_75_toRestore.png';
% name= 'ExtraImages/Lenna_scratch3_toRestore.png';
% name= 'ExtraImages/Lenna_scratch4_toRestore.png';

mask = [181 230 29];       % Yelow mask

%Read the image
I = double(imread(name));

[ni, nj, nC] = size(I);

%TO COMPLETE 1
mask = extractMask(I, mask) ; %mask_img(i,j) == 1 means we have lost information in that pixel
                                      %mask(i,j) == 0 means we have information in that pixel

I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);


%%%Parameters for gradient descent (you do not need for week1)
%param.dt = 5*10^-7;
%param.iterMax = 10^4;
%param.tol = 10^-5;

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

% for each channel 

figure(1)
imshow(I);
title('Before')


% Laplace Equation
 Iinp(:,:,1)=sol_Laplace_Equation_Axb(I_ch1, mask, param);
 Iinp(:,:,2)=sol_Laplace_Equation_Axb(I_ch2, mask, param);
 Iinp(:,:,3)=sol_Laplace_Equation_Axb(I_ch3, mask, param);
    
figure(2)
imshow(Iinp)
title('After');
