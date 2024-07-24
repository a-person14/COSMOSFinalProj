%initial conditions
x0 = 0; %init xpos in meters
y0 = 0; %init ypos
theta0 = 0;
v0 = 0;

%parameters
h = .01; %step size

%init theta array
theta = [theta0];

%Omega array
omega = [pi/2, pi/2, pi/2, pi/2, pi/2, pi/2];

%velocities array
xvelocities = [v0, 1, 2, 1, 1, 1];
yvelocities = [v0, 1, 1, 1, 1, 1];

%positions array
x=[x0];
y=[y0];

%creating the array of theta
for i = 1:6
    theta(i+1) = theta(i) + h*omega(i);
end

%creating the array of x,y positions
for i = 1:6 
    x(i+1) = x(i)+h*xvelocities(i)*cos(theta(i)); 
    y(i+1) = y(i)+h*yvelocities(i)*sin(theta(i));
end

plot(x,y)