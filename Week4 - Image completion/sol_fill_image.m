% Input:
% - im: original image. to be completed with a region from...
% - dbImg: image used to fill image 'im'
% - mask image from im used to denote the hole region
% - mask_extended: dilated mask  

function filled_im = sol_fill_image(im,dbImg,mask,mask_extended)

dbImg = imread(dbImg);
dbImg = double(dbImg)/255;

% Return the most similar region from dbImg to be used to complete im
[im_complet,bb] = find_similar_region(im,dbImg,mask,mask_extended);
%[im_complet,bb] = find_similar_region_difference(im,dbImg,mask,mask_extended);
%[im_complet,bb] = find_similar_region_texture(im,dbImg,mask,mask_extended);


% TODO 1: write the number of states
numPatches = 2;
% end TODO 2

width = bb(3);
height = bb(4);
pos_hole = ceil(bb(1:2));

im_hole = im(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1,:);
mask_im = mask(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1); 

%%
% TODO 2: Fill in CreateGridUGMModel function to create the GM structure
% (4-connected Grid)

edgeStruct = sol_CreateGridUGMModel(height, width, numPatches );
% end TODO 2:

if isempty(edgeStruct)
    error('Function: CreateGridUGMModel has to be implemented');
    filled_im = [];
    return;
end

%% COMPUTING NODE POTENTIALS
%
% TODO 3: Build Unary (node) potencials (aka factors)
%
% Hints:
%  - eq. 2 from paper
%  - see function bwdist to compute Euclidean distance

disp('Computing node potentials...');
k = 0.002;

numNodes = height * width;
Vp = zeros(numNodes,numPatches); %Nodes

dist = bwdist(mask_extended,'euclidean');
distances = (k*dist).^3;
for p=1:numNodes %TODO: For each row in Vp (each row is a node) assign potentials
     %HINT: if i,j belongs to the mask, then it has probabilty 0 of being 
     %      from the canditate image and 1 of being mask
     %HINT: if i,j DOES NOT belongs to the mask, then it has probabilty 1
     %      of being from the canditate image and a probability depending 
     %      of the distance of being mask

    [i,j] = ind2sub( [height, width], p);
        
    if mask_im(i,j) == 1 %if i,j belongs to the mask, then it has probabilty 0 of being from the canditate image and 1 of being mask
        Vp(p,:) = [0 1];
    else %if i,j DOES NOT belongs to the mask, then it has probabilty 1 of being from the canditate image and a probability depending of the distance of being mask
        Vp(p,:) = [1 distances(i,j)];
    end
end
% END TODO 3:

%% COMPUTING EDGE POTENTIALS    
%
% TODO 4: Build pairwise factors
%
% Hint: read last paragraph after eq. 2
%

display('Solving graphical model...');
SSD = sqrt((im_complet(:,:,1) - im_hole(:,:,1)).^2 + (im_complet(:,:,2) - im_hole(:,:,2)).^2 + (im_complet(:,:,3) - im_hole(:,:,3)).^2).^2;
h = 1;
%  HINT: the image is a color image, so every pixel is a vector but
%    diff(p,q) is a scalar (at each pixel), so diff(p,q) should be a 
%    norm. Use L2.

edgePot = zeros(numPatches,numPatches,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
    n1 = edgeStruct.edgeEnds(e,1);
    n2 = edgeStruct.edgeEnds(e,2);    
    
    [i1,j1] = ind2sub( [height, width], n1);
    [i2,j2] = ind2sub( [height, width], n2);
    
    potential = exp( abs(SSD(i1,j1) -  SSD(i2,j2))/ h );
    edgePot(:,:,e) = [1 potential; potential 1];
end
% END TO DO 4

%
% TODO 5: Call inference algorithms to solve the MAP problem
%
% UGM decoding functions:
% - UGM_Decode_LBP
maxOfMarginalsLBPdecode = UGM_Decode_LBP(Vp,edgePot,edgeStruct);
% END TO DO 5
if ~isempty(maxOfMarginalsLBPdecode)
    % Reshape image labels and composite final image
    
    mask_filled = reshape(maxOfMarginalsLBPdecode,[height width]);
    mask(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1) = mask_filled-1;
    
    % Apply Poisson blending
    filled_im = PoissonEditing(im,im_complet,mask,ceil(bb)); 
end

%%
% - UGM_Decode_GraphCut
% maxOfMarginalsGMdecode = UGM_Decode_GraphCut(Vp,edgePot,edgeStruct);
% if ~isempty(maxOfMarginalsGMdecode)
%     % Reshape image labels and composite final image
%     
%     mask_filled = reshape(maxOfMarginalsGMdecode,[height width]);
%     mask(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1) = mask_filled-1;
%     
%     % Apply Poisson blending
%     filled_im = PoissonEditing(im,im_complet,mask,ceil(bb)); 
% end


end

function [im_complet,bb] = find_similar_region(im,complet,mask,mask_gm, scales)

if (~exist('scales','var'))
    [ni_t, nj_t, ~]= size(im);
    [ni_s, nj_s, ~]= size(complet);
    sc = max(ni_t/ni_s, nj_t/nj_s);
    scales=sc*[0.1, 0.2, 0.3, 0.4, 0.5, 0.6 0.7 0.81, 0.90, 1, 1.2, 1.3, 1.4];
end


best_error = inf;
complet_or=complet;
for scl = scales
    complet = imresize(complet_or, scl);
    
    cform = makecform('srgb2lab');
    img_to_inpaint = applycform(im,cform); % the RGB image is transformed to an Lab image

    s = regionprops(mask_gm, 'BoundingBox');
    bb=s.BoundingBox;
    width = s.BoundingBox(3);
    height = s.BoundingBox(4);
    pos_hole = ceil(s.BoundingBox(1:2));

    %computing point that does not penalize 
    %(upper left corner of the bound box of the hole to fill)

    resized_mask = imresize(mask_gm,[size(im,1) size(im,2)]);
    s = regionprops(resized_mask);
    point_no_penalization = ceil(s.Centroid)./[size(im,2) size(im,1)];


    area_to_inpaint = img_to_inpaint(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1,:);
    mask_to_inpaint = mask(pos_hole(1,2):pos_hole(1,2)+height-1,pos_hole(1,1):pos_hole(1,1)+width-1,:);



    m = min([size(complet,1) size(complet,2)]);
    s = floor(([size(complet,1) size(complet,2)]-m)/2);
    complet = complet(s(1)+1:s(1)+m,s(2)+1:s(2)+m,:);
    cform = makecform('srgb2lab');
    dbImg_lab = applycform(complet,cform); % the RGB image is transformed to an Lab image

    size_dbImage = size(complet);
    % The best patch to complete the image is searched
    for x=1:5:size_dbImage(1,2)-width,
        for y=1:5:size_dbImage(1,1)-height,
            crop_area =  dbImg_lab(y:y+height-1,x:x+width-1,:); 
            error = crop_area.*(repmat(mask_to_inpaint==0,[1 1 3]))-area_to_inpaint.*(repmat(mask_to_inpaint==0,[1 1 3]));
            error = sum(error(:).^2)/sum(mask_to_inpaint(:)>0);
            error = error* (abs(x/size_dbImage(1,2)-point_no_penalization(1,1))+abs(y/size_dbImage(1,1)-point_no_penalization(1,2)));
            if error < best_error,
                best_error = error;
                best_crop_area = crop_area;
            end
        end
    end
end    
% Display figure
cform = makecform('lab2srgb');
im_complet =  applycform(best_crop_area,cform); % the Lab image is transformed to an RGB image

end

function [im_complet,bb] = find_similar_region_difference(img_target,img_src,mask,mask_gm, scales)

if (~exist('scales','var'))
    [ni_t, nj_t, ~]= size(img_target);
    [ni_s, nj_s, ~]= size(img_src);
    sc = max(ni_t/ni_s, nj_t/nj_s);
    scales=sc*[0.1, 0.2, 0.3, 0.4, 0.5, 0.6 0.7 0.81, 0.90, 1, 1.2, 1.3, 1.4];
end



%Context is the region between mask and mask extended. i.e. the Xor
%operation between both of them
context = xor(mask, mask_gm);

%crop to the smallest rectangle sorrounding the context area
stats = regionprops(or(mask, mask_gm), 'BoundingBox');
bb = stats.BoundingBox;
ul_corner_src = ceil(bb(1:2));
rect_size = ceil(bb(3:4));

context_crp = context(ul_corner_src(2) : ul_corner_src(2)+ rect_size(2)-1, ...
    ul_corner_src(1):ul_corner_src(1)+ rect_size(1)-1, :);


% the RGB image is transformed to an Lab image
cform_r2l = makecform('srgb2lab');
img_target_lab = applycform(img_target,cform_r2l); 

%crop the image
target_crp = img_target_lab(ul_corner_src(2):ul_corner_src(2)+rect_size(2)-1,...
                            ul_corner_src(1):ul_corner_src(1)+rect_size(1)-1,:);

%Vectorize
vect_target_crp=target_crp(context_crp);

%for each possible scale (see the paper)
min_ssd=inf;
for scl = scales
    img_src_rgb = imresize(img_src, scl);
    
    img_src_lab = applycform(img_src_rgb,cform_r2l); % the RGB image is transformed to a Lab image
    ni_src=size(img_src_lab,1);
    nj_src=size(img_src_lab,2);
    
    %upper-left corner    
    i=1;    
    while i+rect_size(2) <= ni_src
        j=1;
        while j+rect_size(1) <= nj_src
            
            %crop the images
            src_crp    = img_src_lab(i:i+rect_size(2)-1 ,j:j+rect_size(1)-1,:);
            %vectorize
            vect_src_crp=src_crp(context_crp);
            

            %compute SSD
            ssd = sum((vect_target_crp-vect_src_crp).^2);
            if ssd < min_ssd
                region = src_crp;
                min_ssd=ssd
            end
            
            j=j+5;
        end
        i=i+15;
    end
    
    % the Lab image is transformed to an RGB image
    if  (i>1 && j>1)
        cform_l2r = makecform('lab2srgb');
        im_complet =  applycform(region,cform_l2r); 
    end
end
end

function [im_complet,bb] = find_similar_region_texture(img_target,img_src,mask,mask_gm, scales)

if (~exist('scales','var'))
    [ni_t, nj_t, ~]= size(img_target);
    [ni_s, nj_s, ~]= size(img_src);
    sc = max(ni_t/ni_s, nj_t/nj_s);
    scales=sc*[0.1, 0.2, 0.3, 0.4, 0.5, 0.6 0.7 0.81, 0.90, 1, 1.2, 1.3, 1.4];
end



%Context is the region between mask and mask extended. i.e. the Xor
%operation between both of them
context = xor(mask, mask_gm);

%crop to the smallest rectangle sorrounding the context area
stats = regionprops(or(mask, mask_gm), 'BoundingBox');
bb = stats.BoundingBox;
ul_corner_src = ceil(bb(1:2));
rect_size = ceil(bb(3:4));

context_crp = context(ul_corner_src(2) : ul_corner_src(2)+ rect_size(2)-1, ...
    ul_corner_src(1):ul_corner_src(1)+ rect_size(1)-1, :);

%crop the image
target_crp = img_target(ul_corner_src(2):ul_corner_src(2)+rect_size(2)-1,...
                        ul_corner_src(1):ul_corner_src(1)+rect_size(1)-1,:);

%Compute the texture desciptor
tmp =mean(target_crp,3);
[grad_target_j, grad_target_i] = gradient(tmp);
modGrad =sqrt(grad_target_i.^2 + grad_target_j.^2);
target_measure=medfilt2(modGrad, [5 5]);
figure(1)
imshow(target_crp, 'InitialMagnification', 'fit');
hold on
imcontour(context_crp, [1 1])
hold off                    

%Vectorize
vect_target_crp=target_measure(context_crp);

%for each possible scale (see the paper)
min_ssd=inf;
for scl = scales
    img_src_rgb = imresize(img_src, scl);
    
    %Compute the texture desciptor
    tmp =mean(img_src_rgb,3);
    [grad_target_j, grad_target_i] = gradient(tmp);
    modGrad =sqrt(grad_target_i.^2 + grad_target_j.^2);
    src_measure=medfilt2(modGrad, [5 5]);

    ni_src=size(img_src_rgb,1);
    nj_src=size(img_src_rgb,2);
    
    %upper-left corner    
    i=1;    
    
    while i+rect_size(2) <= ni_src
        j=1;
        while j+rect_size(1) <= nj_src
            %figure(2)
            %imshow(img_src_rgb(i:i+rect_size(2)-1 ,j:j+rect_size(1)-1,:), 'InitialMagnification', 'fit');
            %hold on
            %imcontour(context_crp, [1 1])
            %hold off
            
            %crop the images
            src_crp    = src_measure(i:i+rect_size(2)-1 ,j:j+rect_size(1)-1,:);
            %vectorize
            vect_src_crp=src_crp(context_crp);
            

            %compute SSD
            ssd = sum((vect_target_crp-vect_src_crp).^2);
            if ssd < min_ssd
                region = img_src_rgb(i:i+rect_size(2)-1 ,j:j+rect_size(1)-1,:);
                min_ssd=ssd
            figure(3)
            imshow(region, 'InitialMagnification', 'fit');
            hold on
            imcontour(context_crp, [1 1])
            hold off
            pause(0.01);
            end
            
            j=j+5;
        end
        i=i+15;
    end
    
    % the Lab image is transformed to an RGB image   
end
    if  (i>1 && j>1)
    
        im_complet =  region;
    end
end
