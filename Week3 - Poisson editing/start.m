clearvars;
dst = double(imread('lena.png'));
src = double(imread('girl.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('mask_src_eyes.png'));
mask_dst=logical(imread('mask_dst_eyes.png'));

for nC = 1: nChannels
    
    srcI = src(:,:,nC);
    
    %TO DO: COMPLETE the ??
    grad_i = (sol_DiFwd(srcI, param.hi)+sol_DiBwd( srcI, param.hi))./2;%sol_DiBwd( src, param.hi) +  sol_DiFwd( src, param.hi);
    drivingGrad_i = (sol_DiFwd(grad_i, param.hi)+sol_DiBwd( grad_i, param.hi))./2;
    grad_j = (sol_DjFwd(srcI, param.hj)+sol_DjBwd( srcI, param.hj))./2;%sol_DjBwd( src, param.hj) +  sol_DjFwd( src, param.hj);
    drivingGrad_j = (sol_DjFwd(grad_j, param.hj)+sol_DjBwd(grad_j, param.hj))./2;
    
    driving_on_src = drivingGrad_i + drivingGrad_j;%(drivingGrad_i - (2/param.hi).*src)./param.hi + (drivingGrad_j- (2/param.hj).*src)./param.hj ;
    
    driving_on_dst = zeros(size(src(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

%Mouth
%masks to exchange: Mouth
mask_src=logical(imread('mask_src_mouth.png'));
mask_dst=logical(imread('mask_dst_mouth.png'));
for nC = 1: nChannels
    
    srcI = src(:,:,nC);
    
    %TO DO: COMPLETE the ??
    grad_i = (sol_DiBwd( srcI, param.hi)+sol_DiFwd(srcI, param.hi))/2;%sol_DiBwd( src, param.hi) +  sol_DiFwd( src, param.hi);
    drivingGrad_i = (sol_DiBwd(grad_i, param.hi)+sol_DiFwd(grad_i, param.hi))/2;
    grad_j = (sol_DjBwd( srcI, param.hj)+sol_DjFwd(srcI, param.hj))/2;%sol_DjBwd( src, param.hj) +  sol_DjFwd( src, param.hj);
    drivingGrad_j = (sol_DjBwd(grad_j, param.hj)+sol_DjFwd(grad_j, param.hj))/2;
    
    driving_on_src = drivingGrad_i + drivingGrad_j;%(drivingGrad_i - (2/param.hi).*src)./param.hi + (drivingGrad_j- (2/param.hj).*src)./param.hj ;
    
    driving_on_dst = zeros(size(src(:,:,1)));  
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = sol_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)