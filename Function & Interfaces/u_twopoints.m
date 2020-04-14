function r = u_twopoints(x,y)
    [X, Y] = meshgrid(x,y);
    r = min(sqrt((X-1/2).^2+(Y-1/2).^2),sqrt((X+1/2).^2+(Y+1/2).^2));
end