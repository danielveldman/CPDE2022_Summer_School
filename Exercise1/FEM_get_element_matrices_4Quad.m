% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function matrices = FEM_get_element_matrices_4Quad

if exist('FEM_element_matrices_4Quad.mat','file')
   load('FEM_element_matrices_4Quad.mat')
else
    syms x y
    Ne = [(1-x)*(1-y), x*(1-y), x*y, (1-x)*y];
    dxNe = diff(Ne, x);
    dyNe = diff(Ne, y);
    
    matrices.NeNe = double(int(int(Ne'*Ne, x, 0, 1), y, 0, 1));
    matrices.dxNeNe = double(int(int(dxNe'*Ne, x, 0, 1), y, 0, 1));
    matrices.dyNeNe = double(int(int(dyNe'*Ne, x, 0, 1), y, 0, 1));
    matrices.dxNedxNe = double(int(int(dxNe'*dxNe, x, 0, 1), y, 0, 1));
    matrices.dyNedyNe = double(int(int(dyNe'*dyNe, x, 0, 1), y, 0, 1));
    matrices.dxNedyNe = double(int(int(dxNe'*dyNe, x, 0, 1), y, 0, 1));
    
    save('FEM_element_matrices_4Quad.mat', 'matrices');
end