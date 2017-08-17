clearvars;
dst = double(imread('whitewall.jpg'));
src = double(imread('banksy.jpg')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=imread('banksyMaskBefore.jpg');
mask_src_flowers=imread('banksyFlowersMaskBefore.jpg');

% mask_dst=imread('desertMask.jpg');

mask_src = mask_src(:,:,1);
mask_src(mask_src>128) = 256;
mask_src(mask_src<=128) = 0;
mask_src = logical(mask_src);
mask_src(:,1) = 0;

mask_src_flowers = mask_src_flowers(:,:,1);
mask_src_flowers(mask_src_flowers>128) = 256;
mask_src_flowers(mask_src_flowers<=128) = 0;
mask_src_flowers = logical(mask_src_flowers);
mask_src_flowers(:,1) = 0;

mask_dst = padarray(mask_src, [72, 302], 0, 'pre');
mask_dst = padarray(mask_dst, [107,273], 0, 'post');

%Filter Masks
F1=[0 1 0;1 -4 1; 0 1 0];

% monochrome image
srcI = src(:,:,1);
dstI = dst(:,:,1);

grad_dst_i = (sol_DiBwd( dstI, param.hi)+sol_DiFwd(dstI, param.hi))/2;
drivingGrad_dst_i = (sol_DiBwd( grad_dst_i, param.hi)+sol_DiFwd(grad_dst_i, param.hi))/2;
grad_dst_j = (sol_DjBwd( dstI, param.hj)+sol_DjFwd(dstI, param.hj))/2;
drivingGrad_dst_j = (sol_DjBwd( grad_dst_j, param.hj)+sol_DjFwd(grad_dst_j, param.hj))/2;

driving_on_dst_aux = drivingGrad_dst_i + drivingGrad_dst_j;

driving_on_dst = zeros(size(dst(:,:,1)));   
driving_on_dst(mask_dst(:)) = driving_on_dst_aux(mask_dst(:));
    
for nC = 1: nChannels
        
%     dstI = dst(:,:,nC);
    aux_src = src(:,:,nC);
    srcI(mask_src_flowers(:)) = aux_src(mask_src_flowers(:));
    
    img_i = sol_DiFwd(srcI, 1); %i component of the gradient
    img_j = sol_DjFwd(srcI, 1); %j component of the gradient
    
    % div(grad(f))
    grad_v = sol_DiBwd(img_i, 1); %img_ii
    grad_h = sol_DjBwd(img_j, 1); %img_jj
    
    %TO DO: COMPLETE the ??
%     grad_src_i = (sol_DiBwd( srcI, param.hi)+sol_DiFwd(srcI, param.hi))/2;
%     drivingGrad_src_i = (sol_DiBwd( grad_src_i, param.hi)+sol_DiFwd(grad_src_i, param.hi))/2;
%     grad_src_j = (sol_DjBwd( srcI, param.hj)+sol_DjFwd(srcI, param.hj))/2;
%     drivingGrad_src_j = (sol_DjBwd(grad_src_j, param.hj)+sol_DjFwd(grad_src_j, param.hj))/2;
    
    driving_on_src_aux = grad_v + grad_h;

    driving_on_src = zeros(size(dst(:,:,1)));   
    driving_on_src(mask_dst(:)) = driving_on_src_aux(mask_src(:));
    
    param.driving = 1;
    param.driving_dst = 0.8.*driving_on_dst;
    param.driving_src = driving_on_src;

    dstFinalMixing(:,:,nC) = sol_Poisson_Equation_MixingGradients_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dstFinalMixing/256)
