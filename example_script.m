%% Example script 
% Be sure to first run compile_mex_files.m

%% Setup parameters
a = -2;
b = 2;
f = @f_constant;
g = @u_semicircle;
dirichlet = @dirichlet_points_semicircle;
n = 256;

%% Monotone Scheme
tolerance = 10^-6;
um = monotone_solver(n,a,b,dirichlet,g,f,tolerance);

%% Filtered 2nd order centered scheme
p = 1/2;
solver = 'parabolic';
scheme = 'centered';
tolerance = 10^-2;
uf2ndC = filter_solver(n,a,b,dirichlet,g,f,tolerance,p,solver,scheme);

%% Filtered 2nd order upwind scheme
p = 1/2;
solver = 'sweeping';
scheme = '2nd upwind';
tolerance = 10^-3;
uf2ndU = filter_solver(n,a,b,dirichlet,g,f,tolerance,p,solver,scheme);

%% Filtered 2nd order ENO scheme
p = 1/2;
solver = 'sweeping';
scheme = '2nd ENO';
tolerance = 10^-3;
uf2ndENO = filter_solver(n,a,b,dirichlet,g,f,tolerance,p,solver,scheme);
