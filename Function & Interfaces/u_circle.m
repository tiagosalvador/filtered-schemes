function r = u_circle(x,y)
    [X,Y] = meshgrid(x,y);
    r = (X.^2+Y.^2 <= 1).*(0) + (X.^2+Y.^2 > 1).*(sqrt(X.^2+Y.^2)-1);
end