%% test PID code

clear all
close all

%desired (x,y) pos
Setpoints = [[1,1];[2,3]];

%all x,y positions
Px = [0];
Py = [0];

%Velocities (can take this out)
Vx = [1];
Vy = [1];

%step size
h = .1;

%PID x,y gains
Kpx = 2; Kpy = 2; %4,5,4
Kix = .001; Kiy = .001;
Kdx = .01; Kdy = .01;

%init integral+derivatives
xIntegral = 0; yIntegral = 0;
xDerivative = 0; yDerivative = 0;

%xcontroller
for i = 1:length(Setpoints)
    %error initialization
    xError = Setpoints(i,1) - Px(1);
    yError = Setpoints(i,2) - Py(1);
    for k = 1:50 %while loop doesnt work (want to do based on if error is < .0001, move on or do it w time?)
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
        
        plot(Px(k+1),Py(k+1), '-x','MarkerSize',12); drawnow
        hold on
    end 
    % plot(Px,Py)
    % hold on
     Px = [Px(50)];
     Py = [Py(50)];
     Vx = [Vx(50)];
     Vy = [Vy(50)];
     % plot(Px, "green")
     % hold on
     % plot(Py, "blue")
     % xlabel("iterations")
     % ylabel("x or y position")
     % hold on

end
