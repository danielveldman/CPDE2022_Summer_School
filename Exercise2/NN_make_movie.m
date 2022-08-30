% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

clear all
close all
clc

load('Results_1000Iterations')

NT = 101;
I = 64;
N = 2;
time = linspace(0,1,NT);

create_movie = 1; 
if create_movie
    movie_name = 'Learning';  % file name
    vidfile = VideoWriter(movie_name,'MPEG-4');
    vidfile.FrameRate = 10;      % change this number to slow down or speed up the movie
    open(vidfile);
    fig = figure;
    set(fig,'color','w');
end

for ii = 1:NT
    stateii = state1(:,ii);
    stateii = reshape(stateii,2,[]);
    plot(stateii(1,1:I/2), stateii(2,1:I/2), 'bo')
    hold on
    plot(stateii(1,I/2+(1:I/2)), stateii(2,I/2+(1:I/2)), 'ro')
    plot(0, 1,'bx')
    plot(2, 1,'rx')
    hold off
    xlabel 'x_1'
    ylabel 'x_2'
    grid on
    axis equal
    axis([-3, 6, -3, 3])
    title(['t = ', num2str(time(ii)), ' [s]'])
    if create_movie
        frame = getframe(fig);
        writeVideo(vidfile, frame);
    end
    pause(0.1);
end

if create_movie
    close(vidfile);
end