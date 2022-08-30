% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

clear all
close all
clc

N = 2;    % state dimension
I = 64;   % number of data points
NT = 101; % number of time points (layers)
weights = [0.01, 0.01];  % weights [w1, w2]

sigma =@(x) max(x, 0);   % activation function
dsigma =@(x) (x > 0);    % derivative of the activation function

V0 = zeros(2,2,NT-1);    % initial guess for Vk
b0 = zeros(2,NT-1);      % initial guess for bk

time = linspace(0,1,NT); % time grid
dt = time(2) - time(1);

%% generate training data and visualize it
theta = 2*pi*rand(1,I);
r     = [0.4*rand(1,I/2), 0.8+0.4*rand(1,I/2)];
Xin    = [r.*cos(theta); r.*sin(theta)];

figure(1)  % visiualize the input data
plot(Xin(1,1:I/2), Xin(2,1:I/2), 'b.')
hold on
plot(Xin(1,I/2+(1:I/2)), Xin(2,I/2+(1:I/2)), 'r.')
axis equal
xlabel 'x_1'
ylabel 'x_2'
title 'input data'
grid on

Yout = [0*ones(1,I/2), 2*ones(1,I/2);  % final condition in matrix form
    ones(1,I)];

figure(2)  % visiualize the output data
plot(Yout(1,1:I/2), Yout(2,1:I/2), 'b.')
hold on
plot(Yout(1,I/2+(1:I/2)), Yout(2,I/2+(1:I/2)), 'r.')
xlabel 'x_1'
ylabel 'x_2'
title 'output data'
axis equal
grid on

%% a. set up the neural network and loss function
Xin = Xin(:);  % convert the initial data to a column vector
Yout = Yout(:);  % convert the final data to a column vector

state0 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt); % time history

J0 = NN_loss_function(Yout, state0, V0, b0, NT, weights, dt);
Jstart = J0;

%% b. gradient computation and validation
Phi = NN_compute_adjoint(Yout, state0, V0, b0, N, I, NT, dsigma, weights, dt);
[gradV, gradb] = NN_compute_gradients(Phi, state0, V0, b0, N, I, NT, sigma, weights, dt);

% test gradb
db = -gradb;
% db = 0*b0; db(:, end) = rand(1,2);
atest = linspace(0, 0.3, 10);
Jtest = zeros(size(atest));
for ii = 1:length(atest)
    b1 = b0 + atest(ii)*db;
    state1 = NN_compute_state(Xin, V0, b1, N, I, NT, sigma, dt);
    Jtest(ii) = NN_loss_function(Yout, state1, V0, b1, NT, weights, dt);
end
Jappr = Jtest(1) + sum(sum(gradb.*db))*atest;
figure(4)
plot(atest, Jtest, atest, Jappr)
xlabel 'stepsize \beta'
ylabel 'loss function'
title 'J(V_0, b_0 - \beta \nabla_b J)'
grid on
%
% test gradV
dV = -gradV;
atest = linspace(0, 3, 10);
Jtest = zeros(size(atest));
for ii = 1:length(atest)
    V1 = V0 + atest(ii)*dV;
    state1 = NN_compute_state(Xin, V1, b0, N, I, NT, sigma, dt);
    Jtest(ii) = NN_loss_function(Yout, state1, V1, b0, NT, weights, dt);
end
Jappr = Jtest(1) + sum(sum(sum(gradV.*dV)))*atest;
figure(5)
plot(atest, Jtest, atest, Jappr)
xlabel 'stepsize \beta'
ylabel 'loss function'
title 'J(V_0 - \beta \nabla_V J, b_0)'
grid on

%% c. gradient descent
tic
Niters = 100;     % number of iterations (epochs)
beta = 1;         % initial learning rate (learning rate is adapted)
beta_hist = [];
Jhist = J0;
for ii = 1:Niters
    % compute adjoint state and gradients
    
    % choose the learning rate such that the cost function is decreasing
    
    % update values for new iteration
    J0 = J1;
    V0 = V1;
    b0 = b1;
    state0 = state1;
    
    % store history
    Jhist = [Jhist, J0];
    beta_hist = [beta_hist, beta];
end
toc

save(['Results_', num2str(Niters), 'Iterations']) % store results, for example such that you can make a movie later

%% plotting
Xout = state1(:, end);
Xout = reshape(Xout, 2, I);
figure(6)
plot(Xout(1,1:I/2), Xout(2,1:I/2), 'b.')
hold on
plot(Xout(1,I/2+(1:I/2)), Xout(2,I/2+(1:I/2)), 'r.')
xlabel 'x_1'
ylabel 'x_2'
title 'X_{out}'
axis equal
grid on

figure(7)
plot(state1(1:N:I, :).', state1(2:N:I, :).', 'b-')
hold on
plot(state1(I+(1:N:I), :).', state1(I+(2:N:I), :).', 'r-')
plot(state1(1:N:I, end).', state1(2:N:I, end).', 'b.')
plot(state1(I+(1:N:I), end).', state1(I+(2:N:I), end).', 'r.')
xlabel 'x_1'
ylabel 'x_2'
title 'trajectories'
axis equal
x1min = min(state1(1:N:end)); 
x1max = max(state1(1:N:end)); 
x2min = min(state1(2:N:end)); 
x2max = max(state1(2:N:end)); 
axis([1.1*x1min, 1.1*x1max, 1.1*x2min, 1.1*x2max])
grid on

figure(8)
plot(time,state1(1:N:I, :).', 'b-')
hold on
plot(time,state1(I+(2:N:I), :).', 'r-')
ylabel 'x_{i,1}'
xlabel 't (layers)'
grid on

figure(9)
plot(time, state1(2:N:I, :).', 'b-')
hold on
plot(time, state1(I+(2:N:I), :).', 'r-')
ylabel 'x_{i,2}'
xlabel 't (layers)'
grid on

figure(10)
plot(Jhist, 'linewidth', 2)
xlabel 'number of iterations'
ylabel 'loss function'
axis tight
grid on
ylim([0,Inf])

figure(11)
plot(beta_hist)
xlabel 'number of iterations'
ylabel 'learning rate'
axis tight
grid on
ylim([0,Inf])