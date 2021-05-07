%%Define constants
%Wind speed
Vwind_mean = 10;                %m/s
Vwind_standarddev = 3;

%Wind Direction
thetaWind_mean = 0;             %Radians
thetaWind_standarddev = pi/24;

%Particle Diameter
Dp_mean = 65e-6;                %meters
Dp_standarddev = 39e-6;

%Particle Density
rhop_mean = 0.125*1000;         %kg/m^3
rhop_standarddev = 0.005*1000;

%Fluid Properties of air
rhof = 1.225; %kg/m^3

%Parameters of the ship/Blower
htanker = 24.5; %meters
rBlower = 0.762; %meters, 60in/2
vFan = 12.94; %m/s, when cfm = 50000 and d = 60in
vblower = 10; %Fblower/mp; %the initial velocity of the particle, given vFan
%vwind = 10; %meters
%Re = vwind*Dp*rhof/1.68e-5

%Other
Cd = 0.659692;
g = 9.81;

%Number of Particles (iterations)
numParticles = 10;

%Create Output Vectors
xFin = zeros(numParticles,1);
yFin = zeros(numParticles,1);

tFinODE = zeros(numParticles,1);
tFinAnal = zeros(numParticles,1);
err = zeros(numParticles,2);


for i = 1:numParticles
    %%Create Random Parameters
    
    %Random Wind
    %vblower = random('Normal', 20,5);
    vwind = random('Normal', Vwind_mean,Vwind_standarddev);
    thetaWind = random('Normal',thetaWind_mean,thetaWind_standarddev);
    vwindX = vwind*cos(thetaWind);
    vwindY = vwind*sin(thetaWind);
    wind = [vwindX, vwindY];
    
    %Random Diameter
    Dp = abs(random('Normal', Dp_mean, Dp_standarddev));    %made this abs so that it only outputs positive values for the diameter
    Dp = Dp_mean;
    rp = Dp/2;
    A = pi*rp^2;
    Vp = 4/3*pi*rp^3;
    
    %Random Density
    %min and max are roughly 3sd away from mean (min = 0.1, max = 0.14) 
    %source: -3M glass bubbles KS and IM series
    rhop = random('Normal', rhop_mean,rhop_standarddev); 
    mp = rhop*Vp;
    
    FBlower = A*rhof*vFan^2;
    
    %Pack parameters into vector for createState Function
    parameters = zeros(6,1);        %Adjust if you want to add more parameters
    parameters(1) = Cd;
    parameters(2) = A;
    parameters(3) = rhof;
    parameters(4) = mp;
    parameters(5) = g;
    parameters(6) = rBlower;
    parameters(7) = FBlower;
    
    
    %Create Vector of Initial Values
    initX = 0;
    initY = 0;
    initZ = htanker;
    initVx = 0;
    initVy = 0;
    initVz = vblower;
    x0 = [initX; initY; initZ; initVx; initVy; initVz];
    
    
%Adjustment to ode solver to prevent squiggles, increases comp time
    options = odeset('RelTol',1e-5,'AbsTol',1e-7); 
   
%Solving for v of particle v' = .....
    [t,y] = ode23(@(t,x)createState(t,x,parameters,wind),[0,100],x0,options);
    
    xPosition = y(:,1);
    yPosition = y(:,2);
    zPosition = y(:,3);
    xVelocity = y(:,4);
    yVelocity = y(:,5);
    zVelocity = y(:,6);
    
 %Find the index for when particles reach ground.
    [~,groundIndex] = min(abs(zPosition));
 
 %{
 %^^Put a { here to pause state plotting
 
    figure(1)
    subplot(2,1,1)
    plot(t(1:groundIndex), zVelocity(1:groundIndex));
    hold on
    title('z velocity + position')
    xlabel('time')
    ylabel('z-direction velocity')
    subplot(2,1,2)
    plot(t(1:groundIndex), zPosition(1:groundIndex));
    hold on
    xlabel('time')
    ylabel('z-direction position')
    %
    figure(2)
    subplot(2,1,1)
    plot(t(1:groundIndex), xVelocity(1:groundIndex));
    hold on
    title('x velocity + position')
    xlabel('time')
    ylabel('x-direction velocity')
    subplot(2,1,2)
    plot(t(1:groundIndex), xPosition(1:groundIndex));
    hold on
    xlabel('time')
    ylabel('x-direction position')

    figure(3)
    subplot(2,1,1)
    plot(t(1:groundIndex), yVelocity(1:groundIndex));
    hold on
    title('y velocity + position')
    xlabel('time')
    ylabel('y-direction velocity')
    subplot(2,1,2)
    plot(t(1:groundIndex), yPosition(1:groundIndex));
    hold on
    xlabel('time')
    ylabel('y-direction position')

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
%}

%Add to output vector to generate Scatter Plot
xFin(i) = xPosition(groundIndex);
yFin(i) = yPosition(groundIndex);



%
%Analytical Solution
% mp*Vz' - (Cd*A*rhof/2*Vz^2 =  - mp*g);
initVz = 0.5;
a = mp;
b = (Cd*A*rhof)/(2);
c = mp*g;
c1 = 1/sqrt(b*c)*atanh(-initVz*sqrt(b/c));
%c1 = real(c1);
c2 = htanker + a/b*log(cosh(c1*sqrt(b*c)));
%z_analytical = c2 - a/b*log(cosh(sqrt(b*c)/a*(t(1:groundIndex) + a*c1)));
%z_dot = -sqrt(c/b)*tanh(sqrt(b*c)/a*(a*c1 + t(1:groundIndex)));


%Solve for t when z_analytial = 0;
%figure(2)
%plot(t(1:groundIndex),real(z_dot))

%acosh(exp(x) = (apprx) x + 0.693
%t_z0 = a/sqrt(b*c)*acosh(exp(b/a*c2)) - a*c1;
t_z0 = a/sqrt(b*c)*(b/a*c2 + 0.693) - a*c1;
t_z0 = real(t_z0);

tFinAnal(i) = t_z0;
tFinODE(i) = t(groundIndex);


a = Cd*rhof*A/(2*mp);

c1x = 1/(vwindX - initVx);
c2x = initX - 1/a*(vwindX*c1x - log(c1x));
x_analytical = 1/a*(vwindX*(a*t_z0 + c1x) - log(a*t_z0 + c1x)) + c2x;

if vwindY >= 0             %The form of the y-equation changes depending on if vwindY is pos or neg
    c1y =  1/(vwindY - initVy);
    c2y = initY - 1/a*(vwindY*c1y - log(c1y));
    y_analytical = 1/a*(vwindY*(a*t_z0 + c1y) - log(a*t_z0 + c1y)) + c2y;
else
    c1y = 1/(vwindY*sqrt(2))*atanh((vwindY - initVy)/(vwindY*sqrt(2)));
    c2y = initY + 1/a*log(cosh(vwindY*sqrt(2)*c1y));
    y_analytical = c2y + vwindY*t_z0 - 1/a*((vwindY*sqrt(2)*c1y) - 0.693);
    %actually its, c2y + vwindY*t_z0 - 1/a*log(cosh(vwindY*sqrt(2)*c1y) but the
    %log(cosh()) term is too long, so we approximate it.
end






%}
%Checking that the ODE and Analytical solution match up.
ode = [xFin(i), yFin(i)];
ana = [x_analytical , y_analytical];
err(i,:) = (ode - ana)./ode;


end



%xy = num2str(ode)
%xy_anal = num2str(anal)



%{
%Final Scatter Plot
figure(5)
scatter(xFin,yFin)
title('Scatter plot of HGM final location')
xlabel('x distance')
ylabel('y distance')
%}
%Turn off all the holds on the plots
%{
figure(1)
subplot(2,1,1)
hold off
subplot(2,1,2)
hold off

figure(2)
subplot(2,1,1)
hold off
subplot(2,1,2)
hold off

figure(3)
subplot(2,1,1)
hold off
subplot(2,1,2)
hold off

figure(4)
hold off

figure(6)
hold off
%}
%{
err = abs(tFinODE - tFinAnal)./tFinODE;
plot(err,'x')
avg = sum(err)/length(err);
title(['Average Error = ', num2str(avg)])
%}


function dx = createState(t,x,param, wind)
    %Creates the state vector containing info:
    %[Vx,Vy,Vz,ax,ay,az] from the input x = [xpos,ypos,zpos,Vx,Vy,Vz]
    
    Cd = param(1);
    A = param(2);
    rhof = param(3);
    mp = param(4);
    g = param(5);
    rBlower = param(6);
    
    %While Particles are within blower radius, have a contant upward force.
    xPos = x(1);
    yPos = x(2);
    if xPos^2 + yPos^2 <= rBlower^2
        Fblower = param(7);
    else
        Fblower = 0;
    end
    
    Fblower = 0;
    
    vwindX = wind(1);
    vwindY = wind(2);
    
    %Devine Velocity Components
    Vx = x(4);
    Vy = x(5);
    Vz = x(6);
    
    %Define Acceleration Components
    %These are essentially the DEs we're solving
    ax = (-Cd*A*rhof*Vx*abs(Vx)/(2*mp) + Cd*A*rhof*vwindX*Vx + Cd*A*rhof*vwindX*abs(vwindX)/(2*mp));
    ay = (-Cd*A*rhof*Vy*abs(Vy)/(2*mp) + Cd*A*rhof*vwindY*Vy + Cd*A*rhof*vwindY*abs(vwindY)/(2*mp));
    az = (-Cd*A*rhof*Vz*abs(Vz)/(2*mp) - g + Fblower/mp);
    
    %Pack derivatives into vector for ODE solver
    dx = zeros(size(x));
    dx(1) = Vx;
    dx(2) = Vy;
    dx(3) = Vz;
    dx(4) = ax;
    dx(5) = ay;
    dx(6) = az;
end
    
    

