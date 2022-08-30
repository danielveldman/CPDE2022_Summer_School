% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function [gradV, gradb] = NN_compute_gradients(Phi, X, V, b, N, batch_size, NT, sigma, weights, dt, I)

gradV = zeros(size(V));
gradb = zeros(size(b));
s = sigma(X);

for kk = 1:NT-1
    for ii = 1:batch_size
        ind = N*(ii-1)+(1:N);
        gradV(:,:,kk) = gradV(:,:,kk) + dt*Phi(ind,kk)*(s(ind,kk).');
        gradb(:,kk)   = gradb(:,kk)   + dt*Phi(ind,kk);
    end
end
     
gradV = TODO;
gradb = TODO;