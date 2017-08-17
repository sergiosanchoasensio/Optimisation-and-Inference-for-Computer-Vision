function [ u_next ] = multigrid(u,b)
% example: 
%{
clearvars
u_exact = imread('lena.png');
u_exact = double(u_exact(:,:,1))/256;
u = u_exact;
b = zeros(size(u_exact,1)+2);

u = multigrid(u,b);

imshow(u)
imshow(u_exact-u)
%}
    % Smooth
    for i = 1:5
        u = 0.8 * jacobi(u,b) + 0.2 * u;
    end
    
    % Solve the error equation on a coarser grid
    n = size(u, 1) - 1;
    if n>2
        %disp(n)
        r = residual(u, b);
        % Interpolation fine to coarser
        rc = r(1:2:end, 1:2:end); 
        zmat = zeros( round((n+1)/2) );
        ec = multigrid(zmat, rc);
        % Interpolation coarser to fine
        e = zeros(n+1);
        e(1:2:end,1:2:end) = ec;
        e(2:2:end-1,:) = 0.5 * (e(3:2:end,:) + e(1:2:end-2,:));
        e(:,2:2:end-1) = 0.5 * (e(:,3:2:end) + e(:,1:2:end-2));
        u = u + e;
    end
    
    % Post smoothing
    for i = 1:5
        u = 0.8 * jacobi(u, b) + 0.2 * u;
    end
    
    u_next = u;
end 
