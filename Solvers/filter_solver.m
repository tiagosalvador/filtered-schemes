function u = filter_solver(n,a,b,dirichlet,g,f,tolerance,p,solver,scheme)
% Solves the Eikonal equation on [a,b]x[a,b] using a filtered scheme
% on a uniform grid of size n by n.
% Args:
%    dirichlet: function that determines the boundary domain
%    g: boundary condition
%    f: right-hand side of the Eikonal equation
%    tolerance: stopping criteria
%    solver: indicates the type of solver 'parabolic' or 'sweeping'
%    scheme: indicated the accurate scheme to be used. The options are:
%    'centered', '2nd upwind', '2nd upwind', '3rd upwind', '4th upwind', 
%    '2nd ENO', '3rd ENO', '4th ENO'.
% Returns:
%    Solution u.

%% Setup
x = linspace(a,b,n);
y = linspace(a,b,n);
dx = x(2)-x(1);

if strcmp(solver,'parabolic')
    fast = 0;
elseif strcmp(solver,'sweeping')
    fast = 1;
else
    error('Solver must be ''parabolic'' or ''sweeping''')
end

if strcmp(scheme,'centered')
    r = 0; % not used
    order_accurate = 1;
    if strcmp(solver,'sweeping')
        error('The sweeping solver is not available for the filtered centered scheme')
    end
elseif contains(scheme,'upwind')
    r = 0; % not used
    if contains(scheme,'2nd')
        order_accurate = 2;
    elseif contains(scheme,'3rd')
        order_accurate = 3;
    elseif contains(scheme,'4th')
        order_accurate = 4;
    else
        error('Not an available scheme')
    end
elseif contains(scheme,'ENO')
    if contains(scheme,'2nd')
        r = 2;
        order_accurate = 5;
    elseif contains(scheme,'3rd')
        r = 3;
        order_accurate = 6;
    elseif contains(scheme,'4th')
        r = 4;
        order_accurate = 7;
    else
        error('Not an available scheme')
    end
else
    error('Not an available scheme')
end

%% Setting up initial approximation
[dp, Idp, Jdp] = dirichlet(n,x,y);
u0 = initial_solution(n,x,y,g,Idp,Jdp);
[dp, u0] = augment_boundary(dp,u0,g,x,y,n,2);

% u0 = monotone_solver(n,a,b,dirichlet,g,f,10^-6);
% [dp, u0] = augment_boundary(dp,u0,g,x,y,n,2);

if contains(scheme,'ENO')
    C = zeros(n,n);
    NDD = zeros(n,n);
    NDD(:,1) = u0(:,1);
    
    for m=-r:n-r-1
        for k=0:r
            for s=m:m+k-1
                aux = 1;
                for l=m:m+k-1
                    if not(l==s)
                        aux = aux*(-l);
                    end
                end
                C(m+r+1,k+1) = C(m+r+1,k+1) + aux;
            end
            C(m+r+1,k+1) = C(m+r+1,k+1)/factorial(k);
        end
    end
else
    C = zeros(n,n); % not used
    NDD = zeros(n,n); % not used
end

u = u0;
ff = f(x,x);

%% Fixed point iteration

residual = 1;
while residual > tolerance
    uold = u;
    % forces Matlab to deepcopy u into uold. This is necessary given how
    % the C-Mex function loopFilter is implemented
    uold(1) = u(1); 
    u = loopFilter(u,ff,dp,x,y,p,order_accurate,fast,NDD,C,r);
    residual = dx*dx*sum(abs(u(:)-uold(:)));
end


end