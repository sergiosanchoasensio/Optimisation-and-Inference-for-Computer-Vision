clc;
close all;
clearvars;
%clear all;

addpath UGM/KPM
addpath UGM/compiled
addpath UGM/decode
addpath UGM/examples
addpath UGM/infer
addpath UGM/mex
addpath UGM/minConf_2009
addpath UGM/minFunc_2009
addpath UGM/minFunc_2009/logistic
addpath UGM/misc
addpath UGM/sample
addpath UGM/sub
addpath UGM/train

%%% OPEN INPUT IMAGES
disp('Reading the input data...');

%load images
%TODO: Save on the TestSet Folder the image to restore, and read it
file2read = ['TestSet/Image_ToRestore.png'];
I0 = imread(file2read);
I0 = imresize(I0,0.74); %comment if you don't want to downscale the images
I0 = (double(I0)/255);  % image to fill

%load mask
%TODO: Save on the testSet Folder the mask of the image to resotre.
file2read = ['TestSet/maskClean.bmp'];
mask = imread(file2read);
mask = imresize(mask,0.74); %%comment if you don't want to downscale the images
mask = (double(mask)/255);  %mask of the hole to fill

% Creating the mask
% mask = rgb2gray(mask);
mask(mask~=1) = 0;  %remove gray scale, just binary

%Find area where graphical model will be applied
dist = bwdist(mask,'euclidean');
mask_extended = mask;

%try several values
mask_extended(dist<=40 & mask~=1) = 1; % The added area where graphical model will be applied 

% Get Nresults nearest images and use them to fill in missing pixels of query
% image; display results
Nresults = 7;
disp('Finding the most similar images...');
similar_images = find_similar_images(I0,mask,Nresults,'Images/optionalImages.mat');

%%
completed_imgs = cell(Nresults,1);
for j=1:Nresults
    % The image to complete is filled for each of the similar images
    completed_imgs{j,1} = sol_fill_image(I0,similar_images{j,1},mask,mask_extended);
    figure, imshow(completed_imgs{j,1})
     [~,a]=strtok(similar_images{j},'/');
    imwrite(uint8(completed_imgs{j,1}*255),['Results/IMG_0681_',a(2:end)]);
end