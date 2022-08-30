% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function J = OCP_costfunction(T, U, Q, R, time)

J = 0;
for ii = 1:length(time)-1
    dt = time(ii+1) - time(ii);
    J = J + TODO;
end