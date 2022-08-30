% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function [L1, K, M] = FEM_create_matrices_mechanical(Mesh, E, nu, kb, alpha, H)

kb = kb/E/H*(1-nu^2);

%% element matrices
mat = FEM_get_element_matrices_4Quad;

%% Build global matrices

nn = Mesh.nn;
nne = 4;  % number of nodes in element
nne2 = nne^2;
elem_list = Mesh.elem_list;

% L1 = sparse(2*nn, nn);
% K  = sparse(2*nn, 2*nn);

iL1 = zeros(1, 2*nne2*Mesh.ne);
jL1 = zeros(1, 2*nne2*Mesh.ne);
vL1 = zeros(1, 2*nne2*Mesh.ne);

iK = zeros(1, 4*nne2*Mesh.ne);
jK = zeros(1, 4*nne2*Mesh.ne);
vK = zeros(1, 4*nne2*Mesh.ne);

iM = zeros(1, 2*nne2*Mesh.ne);
jM = zeros(1, 2*nne2*Mesh.ne);
vM = zeros(1, 2*nne2*Mesh.ne);

NeNe_v     = mat.NeNe(:);
dxNedxNe_v = mat.dxNedxNe(:);
dyNedyNe_v = mat.dyNedyNe(:);
dxNeNe_v   = mat.dxNeNe(:);
dyNeNe_v   = mat.dyNeNe(:);
dxNedyNe_v = mat.dxNedyNe(:);
dyNedxNe_v = mat.dxNedyNe';
dyNedxNe_v = dyNedxNe_v(:);

index = 0;
for ii = 1:Mesh.ne
    inde = elem_list(ii, 1:4);
    dx   = elem_list(ii, 5);
    dy   = elem_list(ii, 6);
    
%     L1(2*inde-1, inde) = L1(2*inde-1, inde) + mat.dxNeNe*dy;
%     L1(2*inde  , inde) = L1(2*inde  , inde) + mat.dyNeNe*dx;
%     
%     K(2*inde-1, 2*inde-1) = K(2*inde-1, 2*inde-1) + (mat.dxNedxNe*dy/dx + (1-nu)/2*mat.dyNedyNe*dx/dy);
%     K(2*inde-1, 2*inde  ) = K(2*inde-1, 2*inde  ) + (nu*mat.dxNedyNe + (1-nu)/2*mat.dxNedyNe');
%     K(2*inde  , 2*inde-1) = K(2*inde  , 2*inde-1) + (nu*mat.dxNedyNe' + (1-nu)/2*mat.dxNedyNe);
%     K(2*inde  , 2*inde  ) = K(2*inde  , 2*inde  ) + (mat.dyNedyNe*dx/dy + (1-nu)/2*mat.dxNedxNe*dy/dx);
    [ie, je] = ndgrid(inde, inde);
    ie = ie(:);
    je = je(:);
        
    iL1(nne2* 2*index   +(1:nne2)) = 2*ie-1;
    iL1(nne2*(2*index+1)+(1:nne2)) = 2*ie;
    jL1(nne2* 2*index   +(1:nne2)) = je;
    jL1(nne2*(2*index+1)+(1:nne2)) = je;
    vL1(nne2* 2*index   +(1:nne2)) = dxNeNe_v*dy;
    vL1(nne2*(2*index+1)+(1:nne2)) = dyNeNe_v*dx;
    
    iK(nne2* 4*index   +(1:nne2)) = 2*ie-1;
    iK(nne2*(4*index+1)+(1:nne2)) = 2*ie-1;
    iK(nne2*(4*index+2)+(1:nne2)) = 2*ie;
    iK(nne2*(4*index+3)+(1:nne2)) = 2*ie;
    jK(nne2* 4*index   +(1:nne2)) = 2*je-1;
    jK(nne2*(4*index+1)+(1:nne2)) = 2*je;
    jK(nne2*(4*index+2)+(1:nne2)) = 2*je-1;
    jK(nne2*(4*index+3)+(1:nne2)) = 2*je;
    vK(nne2* 4*index   +(1:nne2)) = dxNedxNe_v*dy/dx + (1-nu)/2*dyNedyNe_v*dx/dy + kb*NeNe_v*dx*dy;
    vK(nne2*(4*index+1)+(1:nne2)) = nu*dxNedyNe_v + (1-nu)/2*dyNedxNe_v;
    vK(nne2*(4*index+2)+(1:nne2)) = nu*dyNedxNe_v + (1-nu)/2*dxNedyNe_v;
    vK(nne2*(4*index+3)+(1:nne2)) = dyNedyNe_v*dx/dy + (1-nu)/2*dxNedxNe_v*dy/dx + kb*NeNe_v*dx*dy;
    
    iM(nne2* 2*index   +(1:nne2)) = 2*ie-1;
    iM(nne2*(2*index+1)+(1:nne2)) = 2*ie;
    jM(nne2* 2*index   +(1:nne2)) = 2*je-1;
    jM(nne2*(2*index+1)+(1:nne2)) = 2*je;
    vM(nne2* 2*index   +(1:nne2)) = NeNe_v*dx*dy;
    vM(nne2*(2*index+1)+(1:nne2)) = NeNe_v*dx*dy;
    
    index = index+1;
end

L1 = sparse(iL1, jL1, vL1, 2*nn,   nn);
K  = sparse(iK , jK , vK , 2*nn, 2*nn);
M  = sparse(iM , jM , vM , 2*nn, 2*nn);

L1 = alpha*E*H/(1-nu)*L1;
K  = E*H/(1-nu^2)*K;