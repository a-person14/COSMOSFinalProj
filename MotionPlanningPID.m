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
hold on

% test PID code

%desired (x,y) pos
Setpoints = [[0,-9];[-9.3,-9];[-9.3,-6.4];[-4,-6.4];[-4,-3.8];[9,-3.8];[9,-1];[-1.3,-1];[-1.3,2.6];[-4.6,2.6];[-4.6,4.6];[3.7,4.6];[3.7,6.7];[6,6.7];[6,4.6];[9,4.6];[9,9];[0,9];[0,10]];

%all x,y positions and omega
Px = [0];
Py = [-10];
theta = [pi/2];

%Velocities (can take this out)
Vx = [1];
Vy = [1];
Omega = [0];

%step size
h = .1;

%PID x,y gains
Kpx = 4; Kpy = 2.5; Kptheta = 0.02;%4,5,4 
Kix = .01; Kiy = .1; Kitheta = .05;
Kdx = .01; Kdy = .01; Kdtheta = .04;

%init integral+derivatives
xIntegral = 0; yIntegral = 0; thetaIntegral = 0;
xDerivative = 0; yDerivative = 0; thetaDerivative = 0;

%xcontroller
for i = 1:length(Setpoints)
    %error initialization
    xIntegral = 0; yIntegral = 0; thetaIntegral = 0;
    xDerivative = 0; yDerivative = 0; thetaDerivative = 0;
    xError = Setpoints(i,1) - Px(1);
    yError = Setpoints(i,2) - Py(1);
    k=0;
    thetaError = atan2(yError, xError) - theta(1);
    while norm([xError,yError]) > .08 %while loop doesnt work (want to do based on if error is < .0001, move on or do it w time?)
        k = k+1;
        xIntegral = xIntegral + (Setpoints(i,1) - Px(k)) * h; %integral of the error over time (Setpoints(i,1) - Px(k) here because i wanted to keep prev error in derivation)
        xDerivative = ((Setpoints(i,1) - Px(k)) - xError) / h; %derivative of error in regards to time
        xError = Setpoints(i,1) - Px(k);
    
        %calc new velo
        Vx(k+1) = Kpx*xError + Kix * xIntegral + Kdx*xDerivative; %why is it Kpx*xderivative?
        
        %calc new pos
        Px(k+1) = Px(k) + h*Vx(k)*cos(theta(k));
    
    
        yIntegral = yIntegral + (Setpoints(i,2) - Py(k)) * h;
        yDerivative = ((Setpoints(i,2) - Py(k)) - yError) / h; %(curr err - prev err) / step size
        yError = Setpoints(i,2) - Py(k);
    
        %calc new velo
        Vy(k+1) = Kpy*yError + Kiy * yIntegral + Kdy*yDerivative;
        
        %calc new pos
        Py(k+1) = Py(k) + h*Vy(k)*sin(theta(k));

        %calc theta
        thetaIntegral = thetaIntegral + (atan2(yError, xError) - theta(k)) * h;
        thetaDerivative = ((atan2(yError, xError) - theta(k)) - thetaError) / h;
        thetaError = atan2(yError, xError) - theta(k);
        
        %calc omega
        Omega(k+1) = Kptheta*thetaError;% + Kitheta*thetaIntegral + Kdtheta*thetaDerivative;

        %calc new theta
        theta(k+1) = theta(k)+h*Omega(k);

        plot(Px(k+1),Py(k+1), '.' ,'MarkerSize',5); drawnow
        hold on
    end 
     %plot(Px, "green")
     %hold on
     %plot(Py, "blue")
     %plot(Px,Py)
     %hold on
     Px = [Px(end)];
     Py = [Py(end)];
     Vx = [Vx(end)];
     Vy = [Vy(end)];
     Omega = [Omega(end)];
     theta = [theta(end)];
     % xlabel("iterations")
     % ylabel("x or y position")
     % hold on
     %xlim([-5,10])
     %ylim([-10,-8])

end
