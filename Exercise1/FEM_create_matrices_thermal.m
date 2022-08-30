% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function [CT, KT, HT, VT] = FEM_create_matrices_thermal(Mesh, rho, cp, k, H, Rt, v)

%% element matrices
mat = FEM_get_element_matrices_4Quad;

%% Build global matrices

nn = Mesh.nn;
nne = 4;  % number of nodes in element
nne2 = nne^2;
nnetri = 3;
nnetri2 = nnetri^2;
% nind = Mesh.ne*nne2 + Mesh.netri*nnetri2;
nind = Mesh.ne*nne2;

elem_list = Mesh.elem_list;

% CT = sparse(nn, nn);
% KT = sparse(nn, nn);

iCT = zeros(1, nind);
jCT = zeros(1, nind);
vCT = zeros(1, nind);

vKT = zeros(1, nind);

vVT = zeros(1, nind);

NeNe_v = mat.NeNe(:);
dxNedxNe_v = mat.dxNedxNe(:);
dyNedyNe_v = mat.dyNedyNe(:);
NedxNe_v = mat.dxNeNe';
NedxNe_v = NedxNe_v(:);
NedyNe_v = mat.dyNeNe';
NedyNe_v = NedyNe_v(:);

index = 0;
for ii = 1:Mesh.ne
    inde = elem_list(ii, 1:4);
    dx   = elem_list(ii, 5);
    dy   = elem_list(ii, 6);
    
    %     CT(inde, inde) = CT(inde, inde) + mat.NeNe*dx*dy;
    %
    %     KT(inde, inde) = KT(inde, inde) + mat.dxNedxNe*dy/dx + mat.dyNedyNe*dx/dy;
    
    [ie, je] = meshgrid(inde, inde);
    ie = ie(:);
    je = je(:);
    
    iCT(nne2*index+(1:nne2)) = ie;
    jCT(nne2*index+(1:nne2)) = je;
    vCT(nne2*index+(1:nne2)) = NeNe_v*dx*dy;
    vKT(nne2*index+(1:nne2)) = dxNedxNe_v*dy/dx + dyNedyNe_v*dx/dy;
    vVT(nne2*index+(1:nne2)) = v(1)*NedxNe_v*dy + v(2)*NedyNe_v*dx;
    
    index = index+1;
end

CT = sparse(iCT, jCT, vCT,   nn,   nn);
KT = sparse(iCT, jCT, vKT,   nn,   nn);
VT = sparse(iCT, jCT, vVT,   nn,   nn);

KT = k*H*KT;
HT = 1/Rt*CT;
CT = rho*cp*H*CT;
VT = rho*cp*H*VT;