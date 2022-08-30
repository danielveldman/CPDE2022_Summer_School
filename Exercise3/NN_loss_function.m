% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function J = NN_loss_function(Yout, X, V, b, NT, weights, dt)

J = 0;
for kk = 2:NT-1
    J = J + sum((X(:,kk) - Yout).^2);
end
J = dt*weights(1)*J/2;

J = J + sum((X(:,end) - Yout).^2)/2;

J = J + dt*weights(2)/2*(sum(sum(sum(V.^2))) + sum(sum(b.^2)));