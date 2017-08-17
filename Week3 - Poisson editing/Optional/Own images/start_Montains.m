clearvars;
dst = double(imread('montains.jpg'));
src = double(imread('globo.jpg')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=imread('hotairBalloonMask.jpg');
% mask_dst=imread('desertMask.jpg');

mask_src = mask_src(:,:,1);
mask_src(mask_src>128) = 256;
mask_src(mask_src<=128) = 0;
mask_src = logical(mask_src);
mask_src(:,1) = 0;
mask_src(:,end) = 0;

mask_dst = padarray(mask_src, [358,56], 0, 'pre');
mask_dst = padarray(mask_dst, [0,544], 0, 'post');

%Filter Masks
F1=[0 1 0;1 -4 1; 0 1 0];
    
for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    srcI = src(:,:,nC);
     
    grad_i = (sol_DiBwd( srcI, param.hi)+sol_DiFwd(srcI, param.hi))./2;
    drivingGrad_i = (sol_DiBwd( grad_i, param.hi)+sol_DiFwd(grad_i, param.hi))./2;
    grad_j = (sol_DjBwd( srcI, param.hj)+sol_DjFwd(srcI, param.hj))./2;
    drivingGrad_j = (sol_DjBwd(grad_j, param.hj)+sol_DjFwd(grad_j, param.hj))./2;
    
    driving_on_src = drivingGrad_i + drivingGrad_j;
    
    driving_on_dst = zeros(size(dst(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dstFinal(:,:,nC) = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dstFinal/256)