function r = u_semicircle(x,y)
    [X,Y] = meshgrid(x,y);
    r = (X >= 0 & X.^2+Y.^2 <= 1).*0 + (X >= 0 & X.^2+Y.^2 >= 1).*(sqrt(X.^2+Y.^2)-1) + (X < 0 & abs(Y)<=1).*abs(X) + (X < 0 & Y>1).*sqrt(X.^2+(Y-1).^2) + (X < 0 & Y<-1).*sqrt(X.^2+(Y+1).^2);
end