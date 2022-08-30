% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function X = NN_compute_state(X0, V, b, N, I, NT, sigma, dt)

X = zeros(N*I, NT);
X(:,1) = X0;
for kk = 1:(NT-1)
    for ii = 1:I
        ind = N*(ii-1)+(1:N);
        X(ind,kk+1) = X(ind,kk) + dt*(V(:,:,kk)*sigma(X(ind,kk)) + b(:,kk));
    end
end