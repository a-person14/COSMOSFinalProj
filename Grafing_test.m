clc;
clear;
%close all;
clf

% Define maze boundaries

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
vex = [0 0.5 -0.5 0; 0 -1/sqrt(3) -1/sqrt(3) 0]; % THe vessel's vertecies on the plot
xx= 0
yy= -10
ex = [xx,yy] % the main robot
ex_rotation = 90 % the robot's rotation
ex_rad = deg2rad(ex_rotation) %converts robot's rotation in degrees to radians


sensfldist = 2 %sensorfl offset from main body
sensfltheta = 90  %sensorfl offset angle
sensflthetarad = deg2rad(sensfltheta); % converts offset angle to radians
sensorflx = xx + sensfldist * cos(deg2rad(sensfltheta + ex_rotation)); %robot's sensow front left
sensorfly = yy + sensfldist * sin(deg2rad(sensfltheta + ex_rotation));

sensorfdist = 2 %does sensor front
sensorftheta = 45
sensorfthetarad = deg2rad(sensorftheta)
sensorfx = xx + sensorfdist * cos(deg2rad(sensorftheta + ex_rotation));
sensorfy = yy + sensorfdist * sin(deg2rad(sensorftheta+ ex_rotation));

sensorfrdist = 2
sensorfrtheta = 0
sensorfrthetarad = deg2rad(sensorfrtheta); % the robot's sensor front right
sensorfrx = xx + sensorfrdist * cos(deg2rad(sensorfrtheta + ex_rotation));
sensorfry = yy + sensorfdist * sin(deg2rad(sensorfrtheta + ex_rotation));

front = [xx + cos(ex_rotation), yy + sin(ex_rotation)]

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

plot(xx,yy, front(1), front(2), )
scatter(sensorflx, sensorfly, "filled", "b")
scatter(sensorfx, sensorfy, "filled" , 'g')
scatter(sensorfrx, sensorfry, 'filled', 'o')

wab = 0
% while wab < 1
%     if obstacles{i}(sensorfr(1), sensorfr(2))
%         in_obstacle = true;
%         break;
%     end
%     if in_obstacle == true
%         disp("sensor front right is detecting obstical");
%     end
%     if obstacles{i}(sensorf(1), sensorf(2))
%         in_obstacle = true;
%         break;
%     end
%     if in_obstacle == true
%         disp("sensor front  is detecting obstical");
%     end
%     if obstacles{i}(sensorflx, sensorfly)
%         in_obstacle = true;
%         break;
%     end
%     if in_obstacle == true
%         disp("sensor front left is detecting obstical");
%     end
% 
% end
hold on










% qwe = 1
% the code below makes the bobot move using WASD
% while qwe > 0
%     % Wait for a key press
%     waitforbuttonpress;
% 
%     % Get the current character
%     keyPressed = get(gcf, 'CurrentCharacter');
% 
%     % Display the key pressed
%     disp(['Key pressed: ', keyPressed]);
% 
%     % Break the loop if the 'q' key is pressed
%     if keyPressed == 'w'
%         yy = yy + 1
%         scatter(xx,yy,'filled', 'b')
%     elseif keyPressed == 'a'
%         xx = xx - 1
%         scatter(xx,yy, 'filled', 'b')
% 
%     elseif keyPressed == 'd'
%         xx = xx + 1
%         scatter(xx,yy, "filled", 'b')
%     elseif keyPressed == 's'
%         yy = yy - 1
%         scatter(xx,yy, "filled", 'b')
%     end
% 
% 
% 
% end



% % % % % % % % % % %  THIS IS HOW TO ACCOUNT FOR ROTATION AND STUF
% EHHEHEHAHAHAHAHAHAHAHAHHAHAHAHAHAhAHAHAHHAHAmUAhaHAHAHAHAHHAhAhAHAhAhAHA
% % % % % % % % % % % % Main point coordinates
% % % % % % % % % % % main_x = 0;
% % % % % % % % % % % main_y = 0;   
% % % % % % % % % % % 
% % % % % % % % % % % % Main point rotation (angle in degrees)
% % % % % % % % % % % rotation_angle = 45; % Example rotation angle
% % % % % % % % % % % 
% % % % % % % % % % % % Offset distance and angle from the main point
% % % % % % % % % % % offset_distance = 2;
% % % % % % % % % % % offset_angle = 30; % Angle in degrees
% % % % % % % % % % % 
% % % % % % % % % % % % Convert offset angle to radians
% % % % % % % % % % % offset_angle_rad = deg2rad(offset_angle);
% % % % % % % % % % % 
% % % % % % % % % % % % Calculate subpoint coordinates
% % % % % % % % % % % sub_x = main_x + offset_distance * cos(deg2rad(rotation_angle + offset_angle));
% % % % % % % % % % % sub_y = main_y + offset_distance * sin(deg2rad(rotation_angle + offset_angle));
% % % % % % % % % % % 
% % % % % % % % % % % % Plot main point
% % % % % % % % % % % plot(main_x, main_y, 'ro'); % Red circle for main point
% % % % % % % % % % % 
% % % % % % % % % % % hold on
% % % % % % % % % % % 
% % % % % % % % % % % % Plot subpoint
% % % % % % % % % % % plot(sub_x, sub_y, 'bo'); % Blue circle for subpoint
% % % % % % % % % % % 
% % % % % % % % % % % % Add arrow to show the offset
% % % % % % % % % % % quiver(main_x, main_y, sub_x - main_x, sub_y - main_y, 0, 'b', 'LineWidth', 1.5);
% % % % % % % % % % % 
% % % % % % % % % % % hold off
% % % % % % % % % % % 
% % % % % % % % % % % axis equal









%% make the vessel in matlab, its a triangel
% % % % % Define the vertices of the triangle centered at the origin
% % % % triangle = [0 0.5 -0.5 0; 0 -1/sqrt(3) -1/sqrt(3) 0]; % Triangle vertices
% % % % 
% % % % % Define the rotation angle in degrees
% % % % angle_deg = 30;
% % % % angle_rad = deg2rad(angle_deg); % Convert angle to radians
% % % % 
% % % % % Create the rotation matrix
% % % % R = [cos(angle_rad) -sin(angle_rad); sin(angle_rad) cos(angle_rad)];
% % % % 
% % % % % Rotate and translate the triangle
% % % % rotated_triangle = R * triangle;
% % % % x = rotated_triangle(1, :);
% % % % y = rotated_triangle(2, :);
% % % % 
% % % % % Plot the rotated triangle as a scatter point
% % % % figure;
% % % % fill(x, y, 'b'); % 'b' for blue color, you can choose any color
% % % % axis equal;
% % % % grid on;
% % % % title('Rotated Triangle Scatter Point');
% % % % xlabel('X-axis');
% % % % ylabel('Y-axis');
