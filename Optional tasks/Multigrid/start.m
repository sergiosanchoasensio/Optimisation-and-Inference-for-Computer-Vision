clearvars;
dst = double(imread('gran_canaria.jpg'));
src = double(imread('starwars_tatooine.jpg'));
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

dst1 = dst;

%masks to exchange: Eyes
mask_src=imread('mask_gran_canaria_tatooine.jpg');
mask_src = logical(mask_src>128);
mask_dst=imread('mask_gran_canaria_tatooine.jpg');
mask_dst = logical(mask_dst>128);

se = strel('square',1);
mask_dst = imdilate(mask_dst,se);
mask_src = imdilate(mask_src,se);

[a,b]= find(mask_src(:,:,1)==1); % get the coordinates
x_min = min(a);
x_max = max(a);
y_min = min(b);
y_max = max(b);

bbox = [x_min x_max;y_min y_max];

%For cropping
se = strel('square',1);
mask_dst = imdilate(mask_dst,se);
mask_src = imdilate(mask_src,se);

mask_dst_crop = mask_dst(bbox(1):bbox(3),bbox(2):bbox(4), :);
mask_src_crop = mask_dst(bbox(1):bbox(3),bbox(2):bbox(4), :);
src_crop = src(bbox(1):bbox(3),bbox(2):bbox(4), :);
P = 10;

for nC = 1: nChannels
    
    srcI = src_crop(:,:,nC);
    
    img_i = sol_DiFwd(srcI, 1); %i component of the gradient
    img_j = sol_DjFwd(srcI, 1); %j component of the gradient
    
    % div(grad(f))
    grad_v = sol_DiBwd(img_i, 1); %img_ii
    grad_h = sol_DjBwd(img_j, 1); %img_jj
    
    driving_on_src = grad_v + grad_h;
    driving_on_dst = zeros(size(mask_src_crop(:,:,1)));   
    driving_on_dst = driving_on_src;
    
    param.driving = driving_on_dst;
    
    b = -driving_on_dst;
    %b = zeros(size(mask_src_crop(:,:,1)));
    u = dst(bbox(1):bbox(3),bbox(2):bbox(4), nC);
    
    %{
    for it = 1:50
        u = jacobi(u,b);
    end
    %}
    u = multigrid(u,b);
    
    dst1(bbox(1)+P:bbox(3)-P,bbox(2)+P:bbox(4)-P, nC) = u(P+1:end-P, P+1:end-P);
end

figure;
imshow(dst1/256)