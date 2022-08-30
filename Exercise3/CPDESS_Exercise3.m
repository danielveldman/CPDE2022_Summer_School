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

Vinit = zeros(2,2,NT-1);    % initial guess for Vk
binit = zeros(2,NT-1);      % initial guess for bk

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

%% set up the neural network and loss function
Xin = Xin(:);  % convert the initial data to a column vector
Yout = Yout(:);  % convert the final data to a column vector

state_init = NN_compute_state(Xin, Vinit, binit, N, I, NT, sigma, dt); % time history

Jinit = NN_loss_function(Yout, state_init, Vinit, binit, NT, weights, dt);

Niters = 100;     % number of iterations (epochs)
show_loss = 1;

%% gradient descent with adaptive stepsize
tic
beta = 1;         % initial learning rate (learning rate is adapted)
beta_hist = [];
Jhist = Jinit;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
state0 = state_init;
for ii = 1:Niters
    % compute adjoint state and gradients
    Phi = NN_compute_adjoint(Yout, state0, V0, b0, N, I, NT, dsigma, weights, dt);
    [gradV, gradb] = NN_compute_gradients(Phi, state0, V0, b0, N, I, NT, sigma, weights, dt, I);
    
    % choose the learning rate such that the cost function is decreasing
    J1 = Inf;
    tries = 0;
    beta = beta*4;
    while J1 > J0
        V1 = V0 - beta*gradV;
        b1 = b0 - beta*gradb;
        state1 = NN_compute_state(Xin, V1, b1, N, I, NT, sigma, dt);
        J1 = NN_loss_function(Yout, state1, V1, b1, NT, weights, dt);
        beta = beta / 2;
        tries = tries + 1;
    end
    
    % update values for new iteration
    J0 = J1;
    V0 = V1;
    b0 = b1;
    state0 = state1;
    
    % store history
    if show_loss
        Jhist = [Jhist, J0];
        beta_hist = [beta_hist, beta];
    end
end
toc

if show_loss
    fig = figure(10);
    hold on
    plot(Jhist, 'linewidth', 2)
    xlabel 'number of iterations'
    ylabel 'loss function'
    axis tight
    grid on
    ylim([0,Inf])
    legend('GD (adapt)')
    print(fig, 'gradient_descent_1.jpg', '-djpeg', '-r400');
end

save(['Results_', num2str(Niters), 'Iterations_1'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later

%% 3a. gradient descent with a fixed step size

tic
beta = 0.1;         % initial learning rate
Jhist = Jinit;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
state0 = state_init;
for ii = 1:Niters
    % compute adjoint state and gradients and update
    TODO
    
    % store history
    if show_loss
        state0 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
        J0 = NN_loss_function(Yout, state0, V0, b0, NT, weights, dt);
        Jhist = [Jhist, J0];
    end
end
toc

if show_loss
    figure(10)
    hold on
    plot(Jhist, 'linewidth', 2)
    ylim([0,Inf])
    legend('GD (adapt)', 'GD (fixed)')
    print(fig, 'gradient_descent_2.jpg', '-djpeg', '-r400');
end

state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
V1 = V0;
b1 = b0;
save(['Results_', num2str(Niters), 'Iterations_2'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later

%% stochastic gradient descent
tic
beta = 0.1;         % initial learning rate
Jhist = Jinit;
batch_size = 1;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
state0 = state_init;
for ii = 1:Niters
    
    for jj = 1:I
        batch = randi(I, 1, 1); % select a random data sample
        
        ind = TODO  % compute the indices in the state vector corresponding to the selected sample
        
        state = NN_compute_state(Xin(ind), V0, b0, N, 1, NT, sigma, dt);
        Phi = NN_compute_adjoint(Yout(ind), state, V0, b0, N, 1, NT, dsigma, weights, dt);
        [gradV, gradb] = NN_compute_gradients(Phi, state, V0, b0, N, 1, NT, sigma, weights, dt, I);
        
        V0 = V0 - beta*gradV;
        b0 = b0 - beta*gradb;
    end
    
    % store history
    if show_loss
        state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
        J0 = NN_loss_function(Yout, state1, V0, b0, NT, weights, dt);
        Jhist = [Jhist, J0];
    end
end
toc

if show_loss
    figure(10)
    hold on
    plot(Jhist, 'linewidth', 2)
    ylim([0,Inf])
    legend('GD (adapt)', 'GD (fixed)', 'SGD')
    print(fig, 'gradient_descent_3.jpg', '-djpeg', '-r400');
end

state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
V1 = V0;
b1 = b0;
save(['Results_', num2str(Niters), 'Iterations_3'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later

%% 3c. mini-batch
tic
beta = 0.1;         % initial learning rate
Jhist = Jinit;
batch_size = 4;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
state0 = state_init;
for ii = 1:Niters
    
    for jj = 1:TODO % how many iterations in one epoch
        batch = randi(I, batch_size, 1);
        
        ind = TODO; % select the indices corresponding to the selected data samples
        
        state = NN_compute_state(Xin(ind), V0, b0, N, batch_size, NT, sigma, dt);
        Phi = NN_compute_adjoint(Yout(ind), state, V0, b0, N, batch_size, NT, dsigma, weights, dt);
        [gradV, gradb] = NN_compute_gradients(Phi, state, V0, b0, N, batch_size, NT, sigma, weights, dt, I);
        
        V0 = V0 - beta*gradV;
        b0 = b0 - beta*gradb;
    end
    
    % store history
    if show_loss
        state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
        J0 = NN_loss_function(Yout, state1, V0, b0, NT, weights, dt);
        Jhist = [Jhist, J0];
    end
end
toc

if show_loss
    figure(10)
    hold on
    plot(Jhist, 'linewidth', 2)
    ylim([0,Inf])
    legend('GD (adapt)', 'GD (fixed)', 'SGD', 'minibatch')
    print(fig, 'gradient_descent_4.jpg', '-djpeg', '-r400');
end

state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
V1 = V0;
b1 = b0;
save(['Results_', num2str(Niters), 'Iterations_4'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later

%% 3d. momentum
tic
beta = 0.1;         % initial learning rate
Jhist = Jinit;
batch_size = 1;
gamma = 0.5;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
state0 = state_init;
V1 = V0;
b1 = b0;
for ii = 1:Niters
    
    for jj = 1:I/batch_size
        batch = randi(I, batch_size, 1);
        
        ind = TODO;
        
        state = NN_compute_state(Xin(ind), V1, b1, N, batch_size, NT, sigma, dt);
        Phi = NN_compute_adjoint(Yout(ind), state, V1, b1, N, batch_size, NT, dsigma, weights, dt);
        [gradV, gradb] = NN_compute_gradients(Phi, state, V1, b1, N, batch_size, NT, sigma, weights, dt, I);
        
        V2 = TODO;
        b2 = TODO;
        
        V0 = V1;
        b0 = b1;
        V1 = V2;
        b1 = b2;
    end
    
    % store history
    if show_loss
        state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
        J0 = NN_loss_function(Yout, state1, V0, b0, NT, weights, dt);
        Jhist = [Jhist, J0];
    end
end
toc

if show_loss
    figure(10)
    hold on
    plot(Jhist, 'linewidth', 2)
    ylim([0,Inf])
    legend('GD (adapt)', 'GD (fixed)', 'SGD', 'minibatch', 'momentum')
    print(fig, 'gradient_descent_5.jpg', '-djpeg', '-r400');
end

state1 = NN_compute_state(Xin, V2, b2, N, I, NT, sigma, dt);
V1 = V2;
b1 = b2;
save(['Results_', num2str(Niters), 'Iterations_5'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later

%% 3e. ADAM
tic
beta = 0.1;         % initial learning rate
Jhist = Jinit;
batch_size = 1;
beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;

V0 = Vinit;
b0 = binit;
J0 = Jinit;
mV = 0*V0;
mb = 0*b0;
vV = 0*V0;
vb = 0*b0;
for ii = 1:Niters
    
    for jj = 1:I/batch_size
        batch = randi(I, batch_size, 1);
        
        ind = TODO;
        
        state = NN_compute_state(Xin(ind), V0, b0, N, batch_size, NT, sigma, dt);
        Phi = NN_compute_adjoint(Yout(ind), state, V0, b0, N, batch_size, NT, dsigma, weights, dt);
        [gradV, gradb] = NN_compute_gradients(Phi, state, V0, b0, N, batch_size, NT, sigma, weights, dt, I);
        
        mV = TODO;
        mb = TODO;
        vV = TODO;
        vb = TODO;
        
        kk = (ii-1)*I/batch_size + jj;
        tildemV = TODO;
        tildemb = TODO;
        tildevV = TODO;
        tildevb = TODO;
        
        % update weights and biases
        V0 = TODO;
        b0 = TODO;
    end
    
    %     store history
    if show_loss
        state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
        J0 = NN_loss_function(Yout, state1, V0, b0, NT, weights, dt);
        Jhist = [Jhist, J0];
    end
end
toc

if show_loss
    figure(10)
    hold on
    plot(Jhist, 'linewidth', 2)
    ylim([0,Inf])
    legend('GD (adapt)', 'GD (fixed)', 'SGD', 'minibatch', 'momentum', 'ADAM')
    print(fig, 'gradient_descent_6.jpg', '-djpeg', '-r400');
end

state1 = NN_compute_state(Xin, V0, b0, N, I, NT, sigma, dt);
V1 = V0;
b1 = b0;
save(['Results_', num2str(Niters), 'Iterations_6'], 'state1', 'V1', 'b1') % store results, for example such that you can make a movie later