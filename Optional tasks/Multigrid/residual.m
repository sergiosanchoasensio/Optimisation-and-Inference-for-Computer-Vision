function [ r ] = residual(u,b)
    n = size(u,1) - 1;
    dx = 1.0/n;
    r = zeros(n+1);
    for i = 2:n
        for j = 2:n
            r(i,j) = b(i,j) - ...
                (u(i,j-1) + u(i,j+1) + u(i+1,j) + u(i-1,j) ...
            -4 * u(i,j)) / dx^2;
        end
    end
end
    