clearvars;
I = double(imread('playerMask.jpg'));

I = mean(I,3);

positions = I ~= 0;

I(positions) = 256;

imshow(I)

imwrite(I, 'playerMask.jpg')