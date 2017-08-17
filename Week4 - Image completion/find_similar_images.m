function similar_images = find_similar_images(im,mask,Nresults,image_db)

%% COMPUTING GIST DESCRIPTORS

try
    %All the images of the database have to be preprocesed using the
    %function computeGIST, these images are located in the file
    %descriptorsGIST.m
    %load('descriptorsGIST.mat');
    load(image_db);
catch
    disp(['ERROR: Could not load the .mat file ' GISTmatfile]);
    return;
end
%%
% Color transform structure for sRGB->L*a*b*
cform = makecform('srgb2lab');

% Compute GIST descriptor for the query image
Nfeatures = sum(param.orientationsPerScale)*(param.numberBlocks^2);
if param.color == 1
    imgLab = applycform(im,cform);
    g = zeros(1,3*Nfeatures);
    for j = 1:3
        [des,~] = LMgist(imgLab(:,:,j), '', param);
        g(1,(Nfeatures*(j-1))+1:Nfeatures*j) = des;
    end
else
    [g,~] = LMgist(im, '', param);
end

% Crop both the query image and the mask to the central square portion used
% to compute the GIST descriptor
m = min([size(im,1) size(im,2)]);
s = floor(([size(im,1) size(im,2)]-m)/2);
im = im(s(1)+1:s(1)+m,s(2)+1:s(2)+m,:);
mask = mask(s(1)+1:s(1)+m,s(2)+1:s(2)+m);
mask_sm = imresize(mask,[param.imageSize param.imageSize]);



% GIST weights are the average value of mask pixels within each block
s = floor(linspace(0,param.imageSize,param.numberBlocks+1));
blockWeight = zeros(param.numberBlocks);
for y = 1:param.numberBlocks
    for x = 1:param.numberBlocks
        block = mask_sm(s(y)+1:s(y+1),s(x)+1:s(x+1));
        blockWeight(y,x) = 1 - mean(block(:));
    end
end
blockWeight = repmat(blockWeight(:)',[1 sum(param.orientationsPerScale)]);
if param.color == 1
    blockWeight = repmat(blockWeight,[1 3]);
end

% Calculate weighted Euclidean distance
dist = GIST - repmat(g,[size(GIST,1) 1]);
dist = repmat(blockWeight,[size(GIST,1) 1]).*dist;
dist = sum(dist.^2,2);

% Get Nresults nearest images and use them to fill in missing pixels of query
% image; display results
[~,ind] = sort(dist,'ascend');
best_images = cell(Nresults,1);
for i = 1:Nresults
    % Most similar images are saved
    best_images{i,1} = ['Images' '/' filenames{ind(i)}];
end
similar_images = best_images;
end

