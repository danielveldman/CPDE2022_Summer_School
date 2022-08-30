% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

clear all
close all
clc

%% physical parameters
q     = 1;          % applied heat load [W]
v     = 0.33;       % scanning velocity [m/s]

H     = 0.7e-3;     % substrate thickness [m]
Lslit = 26e-3;      % length of slit (x-direction) [m]
Lfield= 26e-3;      % length of field (x-direction) [m]
Wslit = 5e-3;       % width of slit (y-direction) [m]
Wfield= 33e-3;      % width of field (y-direction) [m]
T     = (Wfield + Wslit)/v; % duration of the simulation [s]

k     = 149;        % substrate thermal conductivity [W/mK]
rho   = 2329;       % substrate massdensity [kg/m3]
cp    = 705;        % substrate heat capacity [J/kgK]
E     = 167e9;      % Young's modulus of substrate [Pa]
nu    = 0.3;        % Poisson's ratio of substrate [-]
kb    = 1.21e12;    % burl stiffness per supported area [N/m3]
alpha = 2.4e-6;     % CTE of substrate [m/mK]
Rt    = 0.9e-3;     % thermal resistance of the burl layer [m2K/W]

NT = 151;             % number of time points
time = linspace(0,T,NT);  % time grid (uniform here)
time2 = time(1:end-1) + diff(time)/2; % centers of the time steps

%% FE model
N1 = 54;
N2 = 32;
N3 = 20;

xgrid = unique([linspace(-4*Lslit, -2*Lslit, N1/4), linspace(-2*Lslit, 2*Lslit , N1), linspace(2*Lslit, 4*Lslit, N1/4)]);
ygrid = unique([linspace(-3*v*T, -2*v*T, N3), linspace(-2*v*T, -1.5*Wslit, N2), linspace(-1.5*Wslit, 1.5*Wslit, N3), linspace(1.5*Wslit, 2*v*T, N2), linspace(2*v*T, 3*v*T, N3)]);

surf(xgrid,ygrid,zeros(length(ygrid), length(xgrid)));
hold on
axis equal
axis tight
view(2)
xlabel 'x [m]'
ylabel '\eta [m]'
title 'used spatial grid'

elemOK =@(x,y) 1;
Mesh = FEM_create_mesh_2D(xgrid, ygrid, elemOK); % create mesh

V = [0, v];    
[CT, KT, HT, VT] = FEM_create_matrices_thermal(Mesh, rho, cp, k, H, Rt, V);   % heat conduction matrices
% [L1, K]          = FEM_create_matrices_mechanical(Mesh, E, nu, kb, alpha, H); % mechanical model

E = CT;
A = -VT -KT -HT;

% disturbance
xr{1} = Lslit/2*[-1 1]; yr{1} = Wslit/2*[-1,1];
Bd  = sum(FEM_create_matrices_mass_part(Mesh, xr{1}, yr{1})).'/Lfield/Wslit;
plot([xr{1}(1), xr{1}(1), xr{1}(2), xr{1}(2), xr{1}(1)], [yr{1}(1), yr{1}(2), yr{1}(2), yr{1}(1), yr{1}(1)], 'b')

% first actuator
xr{2} = Lslit/2*[-1 1]; yr{2} = Wslit/2*[1, 3]; Ar = (xr{2}(2) - xr{2}(1))*(yr{2}(2) - yr{2}(1));
B  = sum(FEM_create_matrices_mass_part(Mesh, xr{2}, yr{2})).' / Ar;
plot([xr{2}(1), xr{2}(1), xr{2}(2), xr{2}(2), xr{2}(1)], [yr{2}(1), yr{2}(2), yr{2}(2), yr{2}(1), yr{2}(1)], 'b')

% second actuator
xr{3} = Lslit/2*[-0.5 0.5]; yr{3} = Wslit/2*[-3, -1]; Ar = (xr{3}(2) - xr{3}(1))*(yr{3}(2) - yr{3}(1));
B(:,2) = sum(FEM_create_matrices_mass_part(Mesh, xr{3}, yr{3})).' / Ar;
plot([xr{3}(1), xr{3}(1), xr{3}(2), xr{3}(2), xr{3}(1)], [yr{3}(1), yr{3}(2), yr{3}(2), yr{3}(1), yr{3}(1)], 'b')

[N, m] = size(B);     % Number of states (N) and number of controls (m)

%% weights for optimal control problem
ind = Mesh.node_list_inv(xr{1}(1) <= xgrid & xgrid <= xr{1}(2), yr{1}(1) <= ygrid & ygrid <= yr{1}(2));
ind = ind(:);
Q = sparse(Mesh.nn, Mesh.nn);
Q(ind, ind) = 1;
R = 0.1*eye(m);

%% 1a. Time discretization by Crank-Nicholson scheme

X0 = zeros(N,1);      % Initial condition (for the state)
U = zeros(m, NT-1);   % Initial guess for the control

T = OCP_compute_temperature(time, X0, E, A, Bd, B, U);

% OPTIONAL: % you can use an decomposition of the matrix E - dt/2*A 
% that needs to be inverted at every time step
% This can speed up computations significantly. 
dt = time(2) - time(1); % uses that the time grid is uniform
[LU.L, LU.U, LU.P, LU.Q, LU.R] = lu(E - dt/2*A);
T = OCP_compute_temperature_octave(time, X0, E, A, Bd, B, U, LU);

xlimits = Lslit*[-1,1]; ylimits = Wfield*[-1.5, 0.5]; 
make_movie = 0; % if you want to create an mp4 movie (1) or not (0). Sorry, you cannot make a movie in  Octave
OCP_display_temperature_field(T, Mesh, time, make_movie, xlimits, ylimits, xr, yr);

%% 1b. Cost function and gradient computation

J = OCP_costfunction(T, U, Q, R, time);

[LUphi.L, LUphi.U, LUphi.P, LUphi.Q, LUphi.R] = lu(E.' - dt/2*A.');
Phi = OCP_compute_adjoint_octave(time, T, E, A, Q, LUphi);
gradU = B.'*Phi + R*U;

figure()
plot(time2, gradU)
xlabel 'time [s]'
ylabel 'gradient of U'

%% 1c. Quadratic approximation

G = OCP_innerproduct(gradU, gradU, time);

dT = OCP_compute_dtemperature_octave(time, E, A, B, gradU, LU);
H = OCP_hessian(dT, gradU, Q, R, time);

beta_opt = G/H;

% verify gradient computation is correct
beta_list = linspace(0,2*beta_opt,20);
for ii = 1:length(beta_list)
    U1 = U - beta_list(ii)*gradU;
    T1 = OCP_compute_temperature_octave(time, X0, E, A, Bd, B, U1, LU);
    J_list(ii) = OCP_costfunction(T1, U1, Q, R, time);
end

figure
plot(beta_list/beta_opt, J_list, 'linewidth', 2)
hold on
plot(beta_list/beta_opt, J - G*beta_list + H/2*beta_list.^2, '--', 'linewidth', 2)
xlabel '\beta / \beta_{opt}'
ylabel 'J'
legend('J(u_0 - \beta \nabla J(u_0))', 'J(u_0) - G\beta + H/2\beta^2')

%% 1d. Gradient descent algorithm

maxiters = 100;
tol      = 1e-5;
for ii = 1:maxiters
    
    % compute gradient (at (U, X))
    
    % compute coefficients in quadratic approximation
    
    % determine the step size
    
%     %%%%%%%%%%%%%%%%%%% CHECK STEPSIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     beta_list = linspace(0,2*beta_opt,20);
%     for ii = 1:length(beta_list)
%         U1 = U - beta_list(ii)*gradU;
%         T1 = OCP_compute_temperature_octave(time, X0, E, A, Bd, B, U1, LU);
%         J_list(ii) = OCP_costfunction(T1, U1, Q, R, time);
%     end
%     
%     figure
%     plot(beta_list/beta_opt, J_list, 'linewidth', 2)
%     hold on
%     plot(beta_list/beta_opt, J - G*beta_list + H/2*beta_list.^2, '--', 'linewidth', 2)
%     xlabel '\beta / \beta_{opt}'
%     ylabel 'J'
%     legend('J(u_0 - \beta \nabla J(u_0))', 'J(u_0) - G\beta + H/2\beta^2')
%     pause
%     %%%%%%%%%%%%%%%%%%% END CHECK STEPSIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % compute new control and cost function value
    U1 = U - beta_opt*gradU;
    T1 = TODO;
    J1 = TODO;
    disp(J1)
    
    % check convergence conditions
    if J - J1 < tol*J
        U = U1;
        T = T1;
        J = J1;
        disp('converged')
        break
    end
    
    % update parameters for next iteration
    U = U1;
    T = T1;
    J = J1;
end

%% plot the obtained results

figure()
plot(time2, U)
xlabel 'time [s]'
ylabel 'u_i(t)'
title 'optimal controls'
legend('u_1(t)', 'u_2(t)')

Bm = 0*B; 
Bm(:,1) = (B(:,1) > 0) / (xr{2}(2) - xr{2}(1)) / (yr{2}(2) - yr{2}(1));
Bm(:,2) = (B(:,2) > 0) / (xr{3}(2) - xr{3}(1)) / (yr{3}(2) - yr{3}(1));
make_movie = 0; % if you want to create an mp4 movie (1) or not (0). Sorry, you cannot make a movie in octave
OCP_display_temperature_and_control(T, Bm*U, Mesh, time, make_movie, xlimits, ylimits, xr, yr);