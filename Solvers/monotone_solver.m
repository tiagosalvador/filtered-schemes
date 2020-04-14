function u = monotone_solver(n,a,b,dirichlet,g,f,tolerance)
    % Solves the Eikonal equation on [a,b]x[a,b] using the monotone scheme
    % on a uniform grid of size n by n.
    % Args:
    %    dirichlet: function that determines the boundary domain
    %    g: boundary condition
    %    f: right-hand side of the Eikonal equation
    %    tolerance: stopping criteria
    % Returns:
    %    Solution u.

%% Setup
x = linspace(a,b,n);
y = linspace(a,b,n);

%% Setting up initial approximation
[dp, Idp, Jdp] = dirichlet(n,x,y);
u = initial_solution(n,x,y,g,Idp,Jdp);

%% Fixed point iteration

ff = f(x,y);
residual = 1;
while residual > tolerance
    uold = repmat(u,1);
    % forces Matlab to deepcopy u into uold. This is necessary given how
    % the C-Mex function loopFilter is implemented
    uold(1) = u(1);
    u = loopMonotone(u,ff,dp,x,y);
    residual = max(max(abs(u-uold)));
end
