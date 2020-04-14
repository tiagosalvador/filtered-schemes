function [dp, Idp, Jdp] = dirichlet_points_two(n,x,y)
    % Determines the boundary indices on the grid when the boundary of the 
    % is the set of points {(1/2,1/2),(-1/2,-1/2)}.
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
points = [1/2 1/2;-1/2 -1/2];
for kk = 1:length(points)
   xcoordinate = points(kk,1); 
   ycoordinate = points(kk,2);
   i = find(y>ycoordinate,1)-1;
   j = find(x>xcoordinate,1)-1;
   dp(i,j) = 1;
   dp(i+1,j) = 1;
   dp(i,j+1) = 1;
   dp(i+1,j+1) = 1;    
end

for i=1:n
    for j=1:n
        if u2(x(j),x(i)) < 0.2
            dp(i,j) = 1;
        end
    end
end

[Idp, Jdp, ~] = find(dp == 1);
end

