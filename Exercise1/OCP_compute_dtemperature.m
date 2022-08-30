% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function dT = OCP_compute_dtemperature(time, E, A, B, dU, LU)

dT = zeros(length(A), length(time));
for ii = 1:length(time)-1
    dt = time(ii+1) - time(ii);
    if nargin >= 6
        % LU factorization available. 
        dT(:,ii+1) = TODO;
    else 
        % LU factorization not available. 
        dT(:,ii+1) = TODO;
    end
end