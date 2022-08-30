% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function out = OCP_innerproduct(u, v, time)

out = 0;
for ii = 1:length(time)-1
    dt = time(ii+1) - time(ii);
    out = out + TODO;
end