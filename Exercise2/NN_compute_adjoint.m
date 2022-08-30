% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function Phi = NN_compute_adjoint(Yout, X, V, b, N, I, NT, dsigma, weights, dt)

Phi = zeros(N*I, NT-1);
Phi(:,end) = TODO;
for kk = NT-2:-1:1
    Vkk = V(:, :, kk+1);
    for ii = 1:I
        ind = N*(ii-1)+(1:N);
        Phi(ind,kk) = Phi(ind,kk+1) + TODO;
    end
end