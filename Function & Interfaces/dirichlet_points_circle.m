function [dp, Idp, Jdp] = dirichlet_points_circle(n,x,y)
    % Determines the boundary indices on the grid when the boundary of the 
    % domain is a circle of radius 1 centered at the origin.
    % Args:
    %    n: size of the grid is nxn
    %    x,y: x and y coordinates of the grid points, respectively
    % Returns:
    %    dp: matrix of size n by n indicating boundary grid points: dp_{ij}
    %    = 1 if the grid point (x_j,y_i) is a boundary grid point; 0
    %    otherwise.
    %    Idp, Jdp: lists of the i and j subscripts of the nonzero indices
    %    of the matrix dp.

dp = zeros(n,n);

for i=1:n
    for j=1:n
        if interface(x(j),y(i))<0
            dp(i,j) = 1;
        end
    end
end

for j=n:-1:1
    s = interface(x(j),y(n)) > 0;
    for i =n:-1:1
        if s ~= (interface(x(j),y(i))>0)
            s = (interface(x(j),y(i))>0);
            dp(i,j) = 1;
            dp(i+1,j) = 1;
            dp(i,j+1) = 1;
            dp(i+1,j+1) = 1;
            dp(i,j-1) = 1;
            dp(i+1,j-1) = 1;
        end
    end
end
for i=n:-1:1
    s = interface(x(n),y(i)) > 0;
    for j = n:-1:1
        if s ~= (interface(x(j),y(i))>0)
            s = (interface(x(j),y(i))>0);
            dp(i,j) = 1;
            dp(i,j+1) = 1;
            dp(i-1,j+1) = 1;
            dp(i+1,j+1) = 1;
            dp(i-1,j) = 1;
            dp(i+1,j) = 1;
        end
    end
end

[Idp, Jdp, ~] = find(dp == 1);

end

function r = interface(x,y)
    [X,Y] = meshgrid(x,y);
    r = sqrt(X.^2+Y.^2)-1;
end