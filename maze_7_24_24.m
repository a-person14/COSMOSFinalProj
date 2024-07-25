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
    if obstacles{i}(1.3, -9.1)
        in_obstacle = true;
        break;
    end
end

if in_obstacle
    disp('(0, 0) is in an obstacle.');
else
    disp('(0, 0) is not in an obstacle.');
end



%%
clc;

h=0.2;
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
esti=10^-2;
x=x0;
y=y0;
t=t0;
w=1;
k_v=0.1;
norm_arr=[goaly-y,goalx-x];
v=k_v*norm(norm_arr);

function in_obstacle = detection (x,y,obstacles)
    in_obstacle = false;
    for i = 1:length(obstacles)
        if obstacles{i}(x, y)
            in_obstacle = true;
            break;
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

function dist = calcdist (x1,y1,x2,y2)
    dist=sqrt((x2-x1)^2+(y2-y1)^2);
end



function path = optipath (posx,posy,pos_arr,t,obstacles,goalx,goaly,h,v,esti)
    path = zeros(3,1);
    pos_x_for=posx + h*v*cos(t); 
    pos_y_for=posy + h*v*sin(t);
    pos_x_left=posx + h*v*cos(t+pi/2); 
    pos_y_left=posy + h*v*sin(t+pi/2);
    pos_x_right=posx + h*v*cos(t-pi/2); 
    pos_y_right=posy + h*v*sin(t-pi/2);
    pos_x_back=posx + h*v*cos(t-pi); 
    pos_y_back=posy + h*v*sin(t-pi);
    is_obs_for=detection(pos_x_for,pos_y_for,obstacles);
    is_obs_left=detection(pos_x_left,pos_y_left,obstacles);
    is_obs_right=detection(pos_x_right,pos_y_right,obstacles);
    for_dist=calcdist(goalx,goaly,pos_x_for,pos_y_for);
    left_dist=calcdist(goalx,goaly,pos_x_left,pos_y_left);
    right_dist=calcdist(goalx,goaly,pos_x_right,pos_y_right);
    pos=[[pos_x_for,pos_y_for,t];[pos_x_left,pos_y_left,t+pi/2];[pos_x_right,pos_y_right,t-pi/2];[pos_x_back,pos_y_back,t-pi]];
    dist=[for_dist,left_dist,right_dist];
    is_obs_arr=[is_obs_for,is_obs_left,is_obs_right];
    min_dist=min(dist);
    max_dist=max(dist);
    min_ind=find(dist==min_dist,1,'last');
    max_ind=find(dist==max_dist,1,'last');
    mid_ind=6-max_ind-min_ind;
    if ((is_obs_arr(min_ind)==0) && (repeatinarr(pos(min_ind,1),pos(min_ind,2),pos_arr,esti)==0))
        path = pos(min_ind,:);
        a=1
    elseif (is_obs_arr(mid_ind)==0) && (repeatinarr(pos(mid_ind,1),pos(mid_ind,2),pos_arr,esti)==0)
        path = pos(mid_ind,:);
        a=2
    elseif (is_obs_arr(max_ind)==0) && (repeatinarr(pos(max_ind,1),pos(max_ind,2),pos_arr(1:2,:),esti)==0)
        path = pos(max_ind,:);
        a=3
    else
        path = pos(4,:);
        a=4
    end
end



i=0;
%while i<5000
while (((abs(pos_arr(1, end)-end_point(1))<esti) && (abs(pos_arr(2, end)-end_point(2))<esti))~=1)
    next_pos=optipath(pos_arr(1,end),pos_arr(2,end),pos_arr,pos_arr(3,end),obstacles,end_point(1),end_point(2),h,v,esti);
    pos_arr(:, end+1) = next_pos;
    i=i+1;
    plot(pos_arr(1,end),pos_arr(2,end),'^','MarkerSize',12); drawnow
    hold on
end