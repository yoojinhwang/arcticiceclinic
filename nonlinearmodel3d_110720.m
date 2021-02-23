%Define constants
Dp = 65e-6; % meters
rp = Dp/2;
A = pi*rp^2;
%rhop = 125; %kg/m^3, 1250 is too large
mp = 4/3*rp^3*pi*rhop;
htanker = 24.5; %meters
rhof = 1.225; %kg/m^3

vFan = 12.94; %m/s, when cfm = 50000 and d = 60in
Fblower = A*rhof*vFan^2;

vblower = 10;   %Fblower/mp; 
                %the initial velocity of the particle, given vFan

%vwind = 10;    %m/s
%Re = vwind*Dp*rhof/1.68e-5;
Cd = 0.659692;
g = 9.81;
time = linspace(0,100,100000);

for i = 1:2
    %vblower = random('Normal', 20,5);
    vwind = random('Normal', 10,3);
    thetaWind = random('Normal',0,pi/4);
    vwindX = vwind*cos(thetaWind);
    vwindY = vwind*sin(thetaWind);
    
    rhop_mean = 0.125*1000; %kg/m^3
    rhop_standarddev = 0.005*1000; %kg/m^3
    %min and max are roughly 3sd away from mean (min = 0.1, max = 0.14) 
    %source: -3M glass bubbles KS and IM series
    rhop = random('Normal', rhop_mean,rhop_standarddev); 
    mp = 4/3*rp^3*pi*rhop;
    
    
    options = odeset('RelTol',1e-5,'AbsTol',1e-7); %adjustment to ode solver to prevent squiggles, increases comp time
    
    %z position
    initZvel = vblower;
    diffeqZ =@(t,zVel)[(-Cd*A*rhof*zVel*abs(zVel)/(2*mp) - g)]; %solving for v of particle v' = .....
    solZ=ode23(diffeqZ,[0,100],initZvel,options);
    plotspace =linspace(0,100,100000);
    zVelocity = deval(solZ,plotspace,1);
    zPosition = time .* zVelocity + htanker;
    
    
 
    for k = 1:length(zPosition)
        if zPosition(k) < 0
            groundIndex = k;
            break
        end
    end
    
    
    
    figure(1)
    subplot(2,1,1)
    plot(time(1:groundIndex), zVelocity(1:groundIndex));
hold on
title('z velocity + position')
xlabel('time')
ylabel('z-direction velocity')
subplot(2,1,2)
plot(time(1:groundIndex), zPosition(1:groundIndex));
hold on
xlabel('time')
ylabel('z-direction position')



    % x position
    initX = 0;
diffeqX = @(t,xVel)[(-Cd*A*rhof*xVel*abs(xVel)/(2*mp) +A*rhof*vwindX*abs(vwindX)/mp)];
solX=ode23(diffeqX,[0,100],initX,options);
plotspace =linspace(0,100,100000);
xVelocity = deval(solX,plotspace,1);
xPosition = time.*xVelocity;

figure(2)
subplot(2,1,1)
plot(time(1:groundIndex), xVelocity(1:groundIndex));
hold on
title('x velocitty + position')
xlabel('time')
ylabel('x-direction velocity')
subplot(2,1,2)
plot(time(1:groundIndex), xPosition(1:groundIndex));
hold on
xlabel('time')
ylabel('x-direction position')


% y position
initY = 0;
diffeqY = @(t,yVel)[(-Cd*A*rhof*yVel*abs(yVel)/(2*mp) +A*rhof*vwindY*abs(vwindY)/mp)];
solY=ode23(diffeqY,[0,100],initY,options);
plotspace =linspace(0,100,100000);
yVelocity = deval(solY,plotspace,1);
yPosition = time.*yVelocity;

figure(3)
subplot(2,1,1)
plot(time(1:groundIndex), yVelocity(1:groundIndex));
hold on
title('y velocity + position')
xlabel('time')
ylabel('y-direction velocity')
subplot(2,1,2)
plot(time(1:groundIndex), yPosition(1:groundIndex));
hold on
xlabel('time')
ylabel('y-direction position')

% x vs z
figure(4)
plot(xPosition(1:groundIndex), yPosition(1:groundIndex))
hold on
xlabel('x position')
ylabel('y position')
title('y vs x')

figure(6)
plot3(xPosition(1:groundIndex),yPosition(1:groundIndex),zPosition(1:groundIndex),'LineWidth',3) 
title('trajectory paths')
xlabel('x distance')
ylabel('y distance')
zlabel('z distance')
hold on
xFin(i) = xPosition(groundIndex);
yFin(i) = yPosition(groundIndex);




end

figure(5)
scatter(xFin,yFin)
title('Scatter plot of HGM final location')
xlabel('x distance')
ylabel('y distance')

%patch([xPosition(1:groundIndex) nan],[yPosition(1:groundIndex) nan],[zPosition(1:groundIndex) nan],[zPosition(1:groundIndex) nan],'EdgeColor','interp','FaceColor','none')
