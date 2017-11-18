clear
clc
close all

% gravitational constant
g = 1;

% start time
tmin = 0;

% start angle (rads)
theta0 = pi / 6;

% frequency of slowest (top) pendulum
f0 = 80;

% number of pendulums
n = 3 * factorial(4);

% frequencies of all pendulums
f = f0: f0 + n - 1;

% periods
T = 1 ./ f;

% lengths
l = (T / (2 * pi)) .^ 2 * g;

% time step
dt = min(T) / 80;

% bounding frame coordinates
ymin = -max(l);
ymax = -min(l) * cos(theta0);
xmin =  ymin * sin(theta0);
xmax = -xmin;

% end time
tmax = 1;

% bounding frame size
dx = xmax - xmin;
dy = ymax - ymin;

% Add margin to frame (15%)
xmax = xmax + 0.15 * dx;
xmin = xmin - 0.15 * dx;
ymax = ymax + 0.15 * dy;
ymin = ymin - 0.15 * dy;

% frame size with margin
dx = xmax - xmin;
dy = ymax - ymin;

% aspect ratio
asp = dx / dy;

% setup colors for each pendulum
colormap('hsv')
c = linspace(0, 1, n);

% video file name
vidobj = VideoWriter('pendulum_wave_12.avi');
vidobj.FrameRate = 60;
vidobj.Quality = 90;
open(vidobj);

fprintf('\nNumber of frames = %d\n', size(tmin: dt: tmax, 2))
fprintf('Aspect ratio = %f\n', asp)

% not sure where I was going with this.  I think I was trying to title each frame with the number of visible pendulum strings (e.g. halfway through the animation, the pendulums split into two distinct string groups, 1/3 and 2/3 of the way through, there are 3 strings, etc.)
nk = 8;
k1 = 1: nk;
k2 = (1: nk)';
k3 = k2 * (1 ./ k1);
k4 = [0; unique(k3(k3 <= 1))]

k = 0;
ip = 0;

% time loop
for t = tmin: dt: tmax
    k = k + 1;
    
    rp = 1000 * (t - tmin) / (tmax - tmin);
    if rp > ip + 1
        ip = floor(rp);
        fprintf('%d / 1000 progress\n', ip);
    end
    
    % angles, x and y coords of each pendulum
    theta = theta0 * cos(2 * pi * t * f);
    x =  l .* sin(theta);
    y = -l .* cos(theta);
    
    % plot
    scatter(x, y, 40, c, 'filled')
    
    whitebg('k')
    axis([xmin xmax ymin ymax])
    pbaspect([asp 1 1])
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    set(gca, 'box', 'off');
    set(gca, 'xcolor', 'k', 'ycolor', 'k')
    set(gcf, 'color', 'k')
    
    %[vmin, imin] = min(abs(k4 - t));
    %title(sprintf('%s', rats(k4(imin))))
    
    fig = gcf;

    % Preserve black bg color.
    fig.InvertHardcopy = 'off';

    % Print to image w/ 247 dpi, i.e. 1080p width.
    cdata = print('-r247', '-RGBImage');
    
    writeVideo(vidobj, cdata);
    
%     if k == 10
%         break
%     end

end

close(vidobj);
