%% test PID code

clear all
close all

%desired (x,y) pos
Setpoints = [[0,0];[1,1]];

%all x,y positions
Px = [0];
Py = [0];

%Velocities (can take this out)
Vx = [1];
Vy = [1];

%step size
h = .01;

%PID x,y gains
Kpx = .6; Kpy = .4;
Kix = .5; Kiy = .5;
Kdx = .4; Kdy = .4;

%init integral+derivatives
xIntegral = 0; yIntegral = 0;
xDerivative = 0; yDerivative = 0;

%xcontroller
for i = 1:length(Setpoints)
    %error initialization
    xError = Setpoints(i,1) - Px(1);
    yError = Setpoints(i,2) - Py(1);
    for k = 1:2000
        xIntegral = xIntegral + (Setpoints(i,1) - Px(k)) * h; %integral of the error over time (Setpoints(i,1) - Px(k) here because i wanted to keep prev error in derivation)
        xDerivative = ((Setpoints(i,1) - Px(k)) - xError) / h; %derivative of error in regards to time
        xError = Setpoints(i,1) - Px(k);
    
        %calc new velo
        Vx(k+1) = Kpx*xError + Kix * xIntegral + Kdx*xDerivative; %why is it Kpx*xderivative?
        
        %calc new pos
        Px(k+1) = Px(k) + h*Vx(k);
    
    
        yIntegral = yIntegral + (Setpoints(i,2) - Py(k)) * h;
        yDerivative = ((Setpoints(i,2) - Py(k)) - yError) / h; %(curr err - prev err) / step size
        yError = Setpoints(i,2) - Py(k);
    
        %calc new velo
        Vy(k+1) = Kpy*yError + Kiy * yIntegral + Kdy*yDerivative;
        
        %calc new pos
        Py(k+1) = Py(k) + h*Vy(k);
        
    end 
    plot(Px, "green")
    hold on
    plot(Py, "blue")
    xlabel("iterations")
    ylabel("x or y position")
    hold on
    Px = [Setpoints(i,1)];
    Py = [Setpoints(i,2)];
end