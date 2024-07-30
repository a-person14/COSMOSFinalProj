clc;
clear;
clf

% Define maze boundaries
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;

% Define obstacles as inequalities in the form x <= f(y)
obstacles = {
    @(x,y) (x >= -8.5 && x <= 7.5) && (y >= -8 && y <= -7.5);
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

% Function to check if a point is in an obstacle
function in_obstacle = is_in_obstacle(x, y, obstacles)
    in_obstacle = false;
    for i = 1:length(obstacles)
        if obstacles{i}(x, y)
            in_obstacle = true;
            break;
        end
    end
end


% Defines a Function To Move Robot Forward
function [x, y] = move_forward(x, y, direction)
    switch direction
        case 0
            y = y + 0.1;
        case 1
            x = x + 0.1;
        case 2 
            y = y - 0.1;
        case 3
            x = x - 0.1;
    end
end


% Defines Function To Detect Walls (to the left) Using Direction in Which
% Robot is Facing
function [front, left] = detect_walls(x, y, direction, obstacles)
    switch direction
        case 0 
            front = is_in_obstacle(x, y+0.1, obstacles);
        case 1
            front = is_in_obstacle(x+0.1, y, obstacles);
        case 2
            front = is_in_obstacle(x, y-0.1, obstacles);
        case 3
            front = is_in_obstacle(x-0.1, y, obstacles);
    end
    
    switch direction
        case 0 
            left = is_in_obstacle(x-0.1, y, obstacles);
        case 1
            left = is_in_obstacle(x, y+0.1, obstacles);
        case 2 
            left = is_in_obstacle(x+0.1, y, obstacles);
        case 3 
            left = is_in_obstacle(x, y-0.1, obstacles);
    end
end


% Robot's Start Position
x = start_point(1);
y = start_point(2);
direction = 0;


% Loop (controls robot movement)
while ~(abs(x - end_point(1)) < 0.1 && abs(y - end_point(2)) < 0.1)
    % Detect walls
    [front, left] = detect_walls(x, y, direction, obstacles);
    
    if left == 0
        % Turn left and move forward
        direction = mod(direction - 1, 4);
        [x, y] = move_forward(x, y, direction);
    elseif front == 0
        % Move forward
        [x, y] = move_forward(x, y, direction);
    else
        % Turn right
        direction = mod(direction + 1, 4);
    end
    
    % Plot robot's position
    plot(x, y, 'bo', 'MarkerSize', marker_size / 2, 'LineWidth', 2);
    pause(0.01);
end


% Displays algorithm is completed
disp('Finished!');
