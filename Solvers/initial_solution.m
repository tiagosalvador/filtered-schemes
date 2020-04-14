function u0 = initial_solution(n,x,y,g,Idp,Jdp)
% Return an initial approximation to the Eikonal equation
% Args:
%    n: grid size
%    x,y: x and y coordinates of the grid points, respectively
%    g: boundary condition
%    Idp,Jdp: subscripts of the boundary indices where the boundary
%    condtion is set.
% Returns:
%    Initial approximate solution u0.

% Initial Approximation
u0 = 10*ones(n,n);

% We impose the Dirichlet boundary
for k = 1:length(Idp)
    u0(Idp(k),Jdp(k)) = g(x(Jdp(k)),y(Idp(k)));
end
end
