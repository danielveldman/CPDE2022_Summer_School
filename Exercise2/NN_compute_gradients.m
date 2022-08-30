% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function [gradV, gradb] = NN_compute_gradients(Phi, X, V, b, N, I, NT, sigma, weights, dt)

gradV = zeros(size(V));
gradb = zeros(size(b));
for kk = 1:NT-1
    for ii = 1:I
        ind = N*(ii-1)+(1:N);
        gradV(:,:,kk) = gradV(:,:,kk) + TODO;
        gradb(:,kk)   = gradb(:,kk)   + TODO;
    end
end