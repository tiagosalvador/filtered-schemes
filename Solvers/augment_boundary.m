function [dp_updated, u0_updated] = augment_boundary(dp,u0,g,x,y,n,p)
% Augments the boundary in order to implement the high-order filtered
% schemes and updates the initial solution accordingly.
% Args:
%    dp: matrix of size n by n indicating the location of the boundary
%    grid points
%    u0: initial solution
%    g: prescribed boundary condition
%    x,y: x and y coordinates of the grid points, respectively
%    p: order of filtered scheme to be used
% Returns:
%    Updated dp, dp_updated, and updated u0, u_updated.

dp_updated = dp;
u0_updated = u0;
for i = 1:n
    for j=1:n
        if dp(i,j)==1
            for k=1:(p-1)
                dp_updated(min(i+k,n),j) = 1;
                dp_updated(max(i-k,1),j) = 1;
                dp_updated(i,max(j-k,1)) = 1;
                dp_updated(i,min(j+k,n)) = 1;
                u0_updated(min(i+k,n),j) = g(x(j),y(min(i+k,n)));
                u0_updated(max(i-k,1),j) = g(x(j),y(max(i-k,1)));
                u0_updated(i,max(j-k,1)) = g(x(max(j-k,1)),y(i));
                u0_updated(i,min(j+k,n)) = g(x(min(j+k,n)),y(i));
            end
        end
    end
end
end