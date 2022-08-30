% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function OCP_display_temperature_field(X, Mesh, time, make_movie, xlimits, ylimits, xr, yr)

duration = 16;
fs = ceil(length(time)/duration);

close all
if make_movie
    vidfile = VideoWriter('movieT.mp4','MPEG-4');
    vidfile.FrameRate = fs;
    open(vidfile);
end

f = figure(1);
set(gca,'nextplot','replacechildren');
for ii = 1:length(time)
    
    surf(Mesh.xgrid, Mesh.ygrid, FEM_shuffle_output_2D(Mesh, X(:, ii)')')
    hold on
    title(['t = ', num2str(time(ii))])
    for kk = 2:length(xr)
        plot3([xr{kk}(1), xr{kk}(1), xr{kk}(2), xr{kk}(2), xr{kk}(1)], ...
            [yr{kk}(1), yr{kk}(2), yr{kk}(2), yr{kk}(1), yr{kk}(1)], 1000*[1 1 1 1 1], 'b-', 'linewidth', 2)
    end
    plot3([xr{1}(1), xr{1}(1), xr{1}(2), xr{1}(2), xr{1}(1)], ...
        [yr{1}(1), yr{1}(2), yr{1}(2), yr{1}(1), yr{1}(1)], 1000*[1 1 1 1 1], 'r--', 'linewidth', 2)
    shading interp
    axis equal
    xlim(xlimits)
    ylim(ylimits)
    xlabel 'x [m]'
    ylabel '\zeta [m]'
    view(2)
    set(gcf,'color','w');
    hold off
    h_cbar = colorbar;
    ylabel(h_cbar, 'temperature increase [K]')
    caxis([min(min(min(X))), max(max(max(X)))])
    
    if make_movie
        frame = getframe(f);
        writeVideo(vidfile,frame);
    else 
        pause(0.1);
    end
end

if make_movie
    close(vidfile);
end