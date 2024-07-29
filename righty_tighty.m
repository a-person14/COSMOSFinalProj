clc;
clear;
%close all;
clf








% the stuff below is all the variables and stuff

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
vex = [0 0.5 -0.5 0; 0 -1/sqrt(3) -1/sqrt(3) 0]; % THe vessel's vertecies on the plot
xx= 0;
yy= -10;

cord = [xx; yy]
bord = 1

ex = [xx,yy] % the main robot
ex_rotation = 0 % the robot's rotation
ex_rad = deg2rad(ex_rotation) %converts robot's rotation in degrees to radians

front = [xx + cos(ex_rotation), yy + sin(ex_rotation)]

sensfldist = 2 %sensorfl offset from main body
sensfltheta = 45  %sensorfl offset angle
sensorflx = xx + sensfldist * cos(sensfltheta + ex_rotation); %robot's sensow front left
sensorfly = yy + sensfldist * sin(sensfltheta + ex_rotation);
% scatter(sensorflx, sensorfly, "filled", "b") does not need to be used??


sensorfdist = 2 %does sensor front
sensorfx = xx + sensorfdist * cos(ex_rotation);
sensorfy = yy + sensorfdist * sin(ex_rotation);
% scatter(sensorfx, sensorfy, "filled" , 'g')


sensorfrdist = 2
sensorfrtheta = -45
sensorfrx = xx + sensorfrdist * cos(sensorfrtheta+ ex_rotation);
sensorfry = yy + sensorfrdist * sin(sensorfrtheta+ ex_rotation);
% scatter(sensorfrx, sensorfry, 'filled', 'o')












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
% in_obstacle = false;
% for i = 1:length(obstacles)
%     if obstacles{i}(-10, 0)
%         in_obstacle = true;
% 
%         break;
%     end
% end
% 
% if in_obstacle %this tells that if an item is in a n obtancle, then it says that it is in an obstace
%     disp('(0, 0) is in an obstacle.');
% else
%     disp('(0, 0) is not in an obstacle.');
% end

hold on






    
qwe = 1
%the code below makes the bobot move using WASD
while qwe > 0
    clc;
    hold on
  scatter(xx,yy, 'filled', 'b')
   front = [xx + cos(ex_rotation), yy + sin(ex_rotation)];
    plot([xx,front(1)], [yy, front(2)], 'LineWidth', 2);
hold on
    end








