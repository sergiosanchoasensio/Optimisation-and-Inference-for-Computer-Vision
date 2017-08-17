function [ u_next ] = jacobi(u,b)
% example: 
%{
u_exact = imread('lena.png');
u_exact = double(u_exact(:,:,1))/256;
u = u_exact;
b = zeros(size(u_exact,1)+2);
for it = 1:100
    u = jacobi(u,b);
end
imshow(u)
imshow(u_exact-u)
%}
    n = size(u,1) - 1;
    dx = 1.0/n;
    u_next = zeros(n+1);
    for i = 2:n
        for j = 2:n
            u_next(i,j) = -dx^2 / 4 * (b(i,j) - ...
                (u(i,j-1) + u(i,j+1) + u(i+1,j) + u(i-1,j)) ...
                / dx^2);
        end
    end
end
                