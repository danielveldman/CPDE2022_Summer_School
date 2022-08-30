% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function Mesh = FEM_create_mesh_2D(xgrid, ygrid, elemOK)

xgrid = sort(unique(xgrid));
ygrid = sort(unique(ygrid));
nx = length(xgrid);
ny = length(ygrid);

nn = nx*ny;           % number of nodes
node_list  = zeros(nn, 2);
node_list_inv = zeros(nx, ny);

ne = (nx-1)*(ny-1);   % number of elements
elem_list = zeros(ne, 6);

kk=0;
mm=0;
for ii = 2:nx
    for jj = 2:ny
        if elemOK(xgrid(ii), ygrid(jj)) && ...
                elemOK(xgrid(ii-1), ygrid(jj)) && ...
                elemOK(xgrid(ii), ygrid(jj-1)) && ...
                elemOK(xgrid(ii-1), ygrid(jj-1))
            
            
            % add nodes
            if ~node_list_inv(ii-1, jj-1)
                kk=kk+1;
                node_list(kk, :) = [ii-1, jj-1];
                node_list_inv(ii-1, jj-1) = kk;
            end
            if ~node_list_inv(ii, jj-1)
                kk=kk+1;
                node_list(kk, :) = [ii, jj-1];
                node_list_inv(ii, jj-1) = kk;
            end
            if ~node_list_inv(ii-1, jj)
                kk=kk+1;
                node_list(kk, :) = [ii-1, jj];
                node_list_inv(ii-1, jj) = kk;
            end
            if ~node_list_inv(ii, jj)
                kk=kk+1;
                node_list(kk, :) = [ii, jj];
                node_list_inv(ii, jj) = kk;
            end
            
            % add element
            nodes = [node_list_inv(ii-1, jj-1), ...
                node_list_inv(ii, jj-1),    ...
                node_list_inv(ii, jj),      ...
                node_list_inv(ii-1, jj) ];
            
            mm=mm+1;
            elem_list(mm, :) = [nodes, ...
                xgrid(ii) - xgrid(ii-1), ...             % dxe (dimensions of the element)
                ygrid(jj) - ygrid(jj-1), ...             % dye
                ];
        end
    end
end
nn = kk;
node_list = node_list(1:nn, :);
ne = mm;
elem_list = elem_list(1:ne, :);


%% store results in a structure called Mesh
Mesh.xgrid         = xgrid;
Mesh.ygrid         = ygrid;
Mesh.nn            = nn;
Mesh.ne            = ne;
Mesh.node_list     = node_list;
Mesh.node_list_inv = node_list_inv;
Mesh.elem_list     = elem_list;
Mesh.netri = 0;
