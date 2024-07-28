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

%%
clc

h=1;
x0=start_point(1);
y0=start_point(2);
t0=0;
goalx=end_point(1);
goaly=end_point(2);
pos_arr=zeros(3,1);
pos_arr(1,1)=x0;
pos_arr(2,1)=y0;
pos_arr(3,1)=t0;
v=1;
t=t0;
esti=10^-2;

function in_obstacle = detection(x, y, obstacles)
    in_obstacle = false;
    if round(x) == 0 && round(y) == 10
        in_obstacle = false;
        return;
    end
    for i = 1:length(obstacles)
        if obstacles{i}(x, y)
            in_obstacle = true;
            return; 
        end
    end
end

function repeat = repeatinarr(x, y, arr, esti)
    repeat = false;
    for i = 1:size(arr, 2) 
        if abs(arr(1, i)-x)<esti && abs(arr(2, i)-y)<esti
            repeat = true;
            break; 
        end
    end
end

function neigh = neighbors (posx,posy,h,v,t)
    pos_x_for=posx + h*v*cos(t); 
    pos_y_for=posy + h*v*sin(t);
    pos_x_left=posx + h*v*cos(t+pi/2); 
    pos_y_left=posy + h*v*sin(t+pi/2);
    pos_x_right=posx + h*v*cos(t-pi/2); 
    pos_y_right=posy + h*v*sin(t-pi/2);
    pos_x_back=posx + h*v*cos(t-pi); 
    pos_y_back=posy + h*v*sin(t-pi);
    neigh = [[pos_x_for,pos_y_for];[pos_x_left,pos_y_left];[pos_x_right,pos_y_right];[pos_x_back,pos_y_back]];
end

function path = bfs(pos_arr, obstacles, start, endpoint, h, v, t, esti)
    rows = 101;
    cols = 101;
    queue = {start};
    visited = false(rows, cols);
    start_ind = round(5*start(1) + 51);
    start_ind_y = round(5*start(2) + 51);
    visited(start_ind, start_ind_y) = true;
    parent = zeros(rows, cols, 2);
    found = false;
    while ~isempty(queue)
        current = queue{1};
        queue(1) = [];
        if isequal(round(current), round(endpoint))
            found = true;
            break;
        end
        next = neighbors(current(1), current(2), h, v, t)
        for i = 1:size(next, 1)
            x_ind = round(5 * next(i, 1) + 51);
            y_ind = round(5 * next(i, 2) + 51);
            if x_ind > 0 && x_ind <= rows && y_ind > 0 && y_ind <= cols
                if ~detection(next(i, 1), next(i, 2), obstacles) && ~visited(x_ind, y_ind)
                    queue{end + 1} = next(i, :);
                    visited(x_ind, y_ind) = true;
                    parent(x_ind, y_ind, 1) = current(1);
                    parent(x_ind, y_ind, 2) = current(2);
                    plot(next(i, 1), next(i, 2), '^', 'MarkerSize', 12); drawnow;
                    hold on;
                end
            end
        end
    end
    
    if found
        path = [];
        node = [endpoint(1),endpoint(2)]; 
        while ~isequal(node, round(start)) 
            path = [node; path]; 
            node(1)
            node(2)
            if node(1) >= -10 && node(1) <= 10 && node(2) >= -10 && node(2) <= 10
               node = squeeze(parent(round(5*node(1)+51), round(5*node(2)+51), :))'; 
            end
        end
        path = [round(start); path];
    else
        path = []; 
    end

end



start_point = [0, -10];
end_point = [0, 10];
path_to_goal = bfs(pos_arr, obstacles, start_point, end_point,h,v,t);


%%
for i=1:size(path_to_goal,1)
    plot(path_to_goal(i,1),path_to_goal(i,2),'*','MarkerSize',12);
end