%Define constants
Dp = 65e-6; % meters
rp = Dp/2;
A = pi*rp^2;
rhop = 1250; %kg/m^3
mp = 4/3*rp^3*pi*rhop;
htanker = 24.5; %meters
rhof = 1.225; %kg/m^3
vFan = 12.94; %m/s, when cfm = 50000 and d = 60in
Fblower = A*rhof*vFan;
vblower = Fblower/mp; %the initial velocity of the particle, given vFan
vwind = 10; %meters
Re = vwind*Dp*rhof/1.68e-5;
Cd = 0; %0.659692;
g = 9.81;
time = linspace(0,40,40000);

for i = 1:1
%     vblower = random('Normal', 15,5);
%     vwind = random('Normal', 8,2);
    
    %z position
    initZvel = vblower;
    diffeqZ =@(t,zVel)[(-Cd*A*rhof*zVel*abs(zVel)/(2*mp) - g)]; %solving for v of particle v' = .....

    solZ=ode45(diffeqZ,[0,40],initZvel);
plotspace =linspace(0,40,40000);
zVelocity = deval(solZ,plotspace,1);
zPosition = time .* zVelocity + htanker;

figure(1)
subplot(2,1,1)
plot(time, zVelocity);
hold on
title('z velocitty + position')
xlabel('time')
ylabel('z-direction velocity')
subplot(2,1,2)
plot(time, zPosition);
hold on
xlabel('time')
ylabel('z-direction position')


% x position
initX = 0;
diffeqX = @(t,xVel)[(-Cd*A*rhof*xVel^2/(2*mp) +A*rhof*vwind^2/mp)];
solX=ode45(diffeqX,[0,40],initX);
plotspace =linspace(0,40,40000);
xVelocity = deval(solX,plotspace,1);
xPosition = time.*xVelocity;

figure(2)
subplot(2,1,1)
plot(time, xVelocity);
hold on
title('x velocitty + position')
xlabel('time')
xlim([0 0.1])
ylabel('x-direction velocity')
subplot(2,1,2)
plot(time, xPosition);
hold on
xlabel('time')
ylabel('x-direction position')


% x vs z
figure(3)
plot(xPosition, zPosition)
hold on
xlabel('x position')
ylabel('z position')
title('z vs x')
ylim([0 30])

end


