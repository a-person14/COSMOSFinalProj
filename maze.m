clc;
clear;
%close all;
clf

% Define maze boundaries
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;

% Define obstacles as inequalities in the form x <= f(y)
obstacles = {
%    @(x,y) y <= x.^2 - 5;                 % Example: parabolic obstacle
%    @(x,y) y >= sin(x) - 3 && y <= 5;     % Example: sine wave obstacle
    @(x,y) (x >= -8.5 && x <= 7.5) && (y >= -8 && y <= -7.5);  % Example: rectangular obstacle
    @(x,y) (x >= 1 && x <= 1.5) && (y >= -10 && y <= -8);
    @(x,y) (x >= 7 && x <= 7.5) && (y >= -7.5 && y <= -5.5);
    @(x,y) (x >= -3 && x <= -2.5) && (y >= -7.5 && y <= -5.5);
    @(x,y) (x >= -3 && x <= 5) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -10 && x <= -8) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -10 && x <= -5) && (y >= -5.5 && y <= -5);
    @(x,y) (x >= -5.5 && x <= -5) && (y >= -5 && y <= -2);
    @(x,y) (x >= -10 && x <= -3) && (y >= 0.5 && y <= 1);
    @(x,y) (x >= -3 && x <= -2.5) && (y >= -2 && y <= 1);
    @(x,y) (x >= -3 && x <= 7.5) && (y >= -2.5 && y <= -2);
    @(x,y) (x >= 0 && x <= 7.5) && (y >= 0.5 && y <= 1);
    @(x,y) (x >= 0 && x <= 0.5) && (y >= 1 && y <= 4);
    @(x,y) (x >= -3 && x <= 10) && (y >= 3.5 && y <= 4);
    @(x,y) (x >= -10 && x <= -6) && (y >= 3.5 && y <= 4);
    @(x,y) (x >= -7.5 && x <= 7.5) && (y >= 7.5 && y <= 8);
    @(x,y) (x >= -1.5 && x <= -1) && (y >= 8 && y <= 10);
    @(x,y) (x >= 7 && x <= 7.5) && (y >= 6 && y <= 7.5);
    @(x,y) (x >= 4.5 && x <= 5) && (y >= 4 && y <= 6);
    @(x,y) (x >= -3 && x <= 1) && (y >= 5.5 && y <= 6);
    @(x,y) (x >= -8 && x <= -7.5) && (y >= 6 && y <= 8);
    @(x,y) (x >= -8 && x <= -2) && (y >= 5.5 && y <= 6);
    @(x,y) x <= -10; 
    @(x,y) x >= 10; 
    @(x,y) y <= -10; 
    @(x,y) y >= 10; 
    
};

% Define start and end points
start_point = [0, -10];
end_point = [0, 10];

% Visualization parameters
resolution = 0.1;  % Resolution for plotting obstacles
marker_size = 10;

% Plot maze boundaries
figure(1);
hold on;
plot([-1, xmin, xmin, -1], [ymax, ymax, ymin, ymin], 'k-', 'LineWidth', 2);
plot([1, xmax, xmax, 1], [ymin, ymin, ymax, ymax], 'k-', 'LineWidth', 2);
%plot([xmin, -1, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k-', 'LineWidth', 2);

% Plot obstacles
x_range = xmin:resolution:xmax;
y_range = ymin:resolution:ymax;
for i = 1:length(obstacles)
    for j = 1:length(x_range)
        for k = 1:length(y_range)
            if obstacles{i}(x_range(j), y_range(k)) && x_range(j)>-10 && x_range(j)<10 && y_range(k)>-10 && y_range(k)<10
                plot(x_range(j), y_range(k), 'r.', 'MarkerSize', 10);
            end
        end
    end
end

% Plot start and end points
plot(start_point(1), start_point(2), 'go', 'MarkerSize', marker_size, 'LineWidth', 2);  % Green circle for start
plot(end_point(1), end_point(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);  % Red circle for end

% Adjust axes and grid
axis equal;
grid on;
xlim([xmin, xmax]);
ylim([ymin, ymax]);
title('COSMOS 24 Cluster 11 Maze');

%% Example of how to check if a point is in an obstacle
% Check if (0, 0) is in an obstacle
in_obstacle = false;
for i = 1:length(obstacles)
    if obstacles{i}(-10, 0)
        in_obstacle = true;
        break;
    end
end

if in_obstacle
    disp('(0, 0) is in an obstacle.');
else
    disp('(0, 0) is not in an obstacle.');
end