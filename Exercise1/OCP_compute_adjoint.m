% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function Phi = OCP_compute_adjoint(time, T, E, A, Q, LU)

N = size(A,1); % dimension of the state space

Phi = zeros(N, length(time)-1);
dt = time(end) - time(end-1);
Phi(:,end) = TODO;
for ii = length(time)-1:-1:2
    dt = time(ii) - time(ii-1);
    if nargin >= 6
        % LU factorization available. 
        Phi(:,ii-1) = TODO;
    else
        % LU factorization not available. 
        Phi(:,ii-1) = TODO;
    end
end