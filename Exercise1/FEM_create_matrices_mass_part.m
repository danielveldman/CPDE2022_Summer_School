% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function E = FEM_create_matrices_mass_part(Mesh, xr, yr)

%% Build global matrices

nn = Mesh.nn;
elem_list = Mesh.elem_list;
node_list = Mesh.node_list;

E  = sparse(nn, nn);

% syms x y
% Ne = [(1-x)*(1-y), x*(1-y), x*y, (1-x)*y];
% intNeNe = int(int(transpose(Ne)*Ne,x),y);

intNeNe =@(x, y) ...
[          ((x - 1)^3*(y - 1)^3)/9,    -(x^2*(2*x - 3)*(y - 1)^3)/18, (x^2*y^2*(2*x - 3)*(2*y - 3))/36,    -(y^2*(2*y - 3)*(x - 1)^3)/18
     -(x^2*(2*x - 3)*(y - 1)^3)/18,        (x^3*y*(y^2 - 3*y + 3))/9,          -(x^3*y^2*(2*y - 3))/18, (x^2*y^2*(2*x - 3)*(2*y - 3))/36
  (x^2*y^2*(2*x - 3)*(2*y - 3))/36,          -(x^3*y^2*(2*y - 3))/18,                      (x^3*y^3)/9,          -(x^2*y^3*(2*x - 3))/18
     -(y^2*(2*y - 3)*(x - 1)^3)/18, (x^2*y^2*(2*x - 3)*(2*y - 3))/36,          -(x^2*y^3*(2*x - 3))/18,  y^3*(x^3/9 - x^2/3 + x/3 - 1/9)]; 

for ii = 1:Mesh.ne
    inde = elem_list(ii, 1:4);
    dx   = elem_list(ii, 5);
    dy   = elem_list(ii, 6);
    
    xe = Mesh.xgrid(node_list(inde, 1));
    xe = [min(xe), max(xe)];
    ye = Mesh.ygrid(node_list(inde, 2));
    ye = [min(ye), max(ye)];
    
    if xe(2) >= xr(1) && xe(1) <= xr(2) && ye(2) >= yr(1) && ye(1) <= yr(2)
        
        xi1  = max(0, min(1, (xr(1) - xe(1)) / dx));
        xi2  = max(0, min(1, (xr(2) - xe(1)) / dx));
        
        et1  = max(0, min(1, (yr(1) - ye(1)) / dy));
        et2  = max(0, min(1, (yr(2) - ye(1)) / dy));
        
        NeNe     = intNeNe(xi2, et2)     - intNeNe(xi2, et1)     - intNeNe(xi1, et2)     + intNeNe(xi1, et1);
    
        E(inde, inde) = E(inde, inde) + NeNe*dx*dy;
    end
end