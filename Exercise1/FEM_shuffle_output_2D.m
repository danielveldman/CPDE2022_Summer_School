% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function out = FEM_shuffle_output_2D(Mesh, values)

nx = length(Mesh.xgrid);
ny = length(Mesh.ygrid);
nt = size(values, 1);
out = zeros(nx, ny, nt);
for ii = 1:nx
    for jj = 1:ny
        ind = Mesh.node_list_inv(ii, jj);
        if ~ind
            out(ii, jj, :) = NaN;
        else
            out(ii, jj, :) = values(:, ind);
        end
    end
end