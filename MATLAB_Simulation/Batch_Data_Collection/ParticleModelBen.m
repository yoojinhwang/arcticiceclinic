function OutData = ParticleModelBen(numParticles,TimeDep,modDenom,stdevWindChange)
    %numParticles: Number of iterations/particles
    %TimeDep: Triggers time dependence (1 on, 0 off)
    %modDenom: Only relevant in TD models, frequency of wind change in
    %seconds.  Must be an integer
    %stdevWindChange: Only relevant in TD models, standard deviation of
    %random wind speed changes (gaussian distribution)
    
    %OutData: A table of output data.  Can be configured to fit your data
    %needs
    
    startTime = tic; %start timer
    
    %Time Dependend Wind Variables
    global windChangeCounter wind WindSpeedTracker

    %%Define constants
    %Wind speed
    Vwind_mean = 10;                %m/s
    Vwind_standarddev = 3;

    %Wind Direction
    thetaWind_mean = 0;             %Radians
    thetaWind_standarddev = pi/24;

    %Particle Diameter
    k = 4.1; %Parameters tuned for 3M K1 HGM
    th = 17e-6;
    Dp_max = 420e-6; %maximum particle diameter
    Dp_min = 10e-6; %minimum particle diameter

    %Particle Density
    rhop_mean = 0.125*1000;         %kg/m^3
    rhop_standarddev = 0.005*1000;

    %Fluid Properties of air
    rhof = 1.341; %kg/m^3 calculated for -10 degrees C
    
    %Parameters of the ship/Blower
    htanker = 24.5; %meters
    rBlower = 0.762; %meters, 60in/2
    vFan = 12.94; %m/s, when cfm = 50000 and d = 60in
    vblower = 10; %Fblower/mp; %the initial velocity of the particle, given vFan
    %vwind = 10; %meters
    %Re = vwind*Dp*rhof/1.68e-5;

    %Other
    Cd = 0.659692; %2.6; FURTHER WORK: Vary CD with reynolds number (pretty well behaved eqn for spheres)
    g = 9.81;
    failedCount=0;

    %Create Output Vectors
    xFin = zeros(numParticles,1);
    yFin = zeros(numParticles,1);
    zFin = zeros(numParticles,1);
    tFin = zeros(numParticles,1);
    Ds = zeros(numParticles,1);
    rhos = zeros(numParticles,1);
    mps = zeros(numParticles,1);
    AvgWindX = zeros(numParticles,1);
    AvgWindY = zeros(numParticles,1);
    AvgTheta = zeros(numParticles,1);
    AvgMag = zeros(numParticles,1);
    %Timing for debug
    times = zeros(numParticles,1);


    for i = 1:numParticles
        %%Create Random Parameters
        runitime=tic;
        %Random Wind
        %vblower = random('Normal', 20,5);
        vwind = random('Normal', Vwind_mean,Vwind_standarddev);
        thetaWind = random('Normal',thetaWind_mean,thetaWind_standarddev);
        wind(1) = vwind*cos(thetaWind); %vwindX
        wind(2) = vwind*sin(thetaWind); %vwindY
        WindSpeedTracker = wind; %resets WindSpeedTracker variable and loads first wind speed


        %Random Diameter
        %Dp = random('Normal',Dp_mean,Dp_standarddev); %Use for standard gaussian distribution (you'll have to make some changes to the if statements if you still want to incorporate min/max cutoffs)
        DpOK=0;
        while DpOK ~=1    
            Dp =random('Gamma', k, th);
            if Dp<Dp_max && Dp>Dp_min
                DpOK=1;
            end
        end
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
        parameters(8) = modDenom;
        parameters(9) = stdevWindChange;

        %Create Vector of Initial Values
        initX = 0;
        initY = 0;
        initZ = htanker;
        initVx = 0;
        initVy = 0;
        initVz = vblower;
        x0 = [initX; initY; initZ; initVx; initVy; initVz];
        
        %Set how many seconds the ODE's are solved to.  Increase if particles aren't hitting the ground
        %Based on testing, critical ODE time limit is almost entirely determined by particle diameter, with smaller particles requiring more time
        if Dp>50e-6
            ODETimeLim = 200; %should bo 100, 200
        else
            ODETimeLim = 300;
        end

        %Adjustment to ode solver to prevent squiggles, increases comp time
        options = odeset('RelTol',1e-5,'AbsTol',1e-7); 

        %Solving for v of particle v' = .....
        windChangeCounter = 0; %Ensures wind changes just once per time interval
        [t,y] = ode23(@(t,x)createState(t,x,parameters,TimeDep),[0,ODETimeLim],x0,options);
        xPosition = y(:,1);
        yPosition = y(:,2);
        zPosition = y(:,3);
        xVelocity = y(:,4);
        yVelocity = y(:,5);
        zVelocity = y(:,6);

        %Find the index for when particles reach ground.
        [minZ,groundIndex] = min(abs(zPosition));
        %Check to see if they actually hit the ground or just ran out of
        %time
        %FURTHER WORK: Find average z velocity of particles in order to
        %determine how many meters of error are possible for a given
        %failure height (eg, asuume vx = 10 m/s, find time (t) it takes
        %for a particle to fall from minZ to the ground, t*vx = final
        %position error).  Use said limit to impelemet a acceptable error
        %threshold instead of the arbritrary choice I made
        if minZ>0.1 %This is an estimation, not an actual count
            fprintf('Error: Particle may not have hit ground, Zpos = %.2f \n',minZ)
            failedCount = failedCount+1;
        end
        
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
        
        %Output Data Collection
        xFin(i) = xPosition(groundIndex);
        yFin(i) = yPosition(groundIndex);
        zFin(i) = zPosition(groundIndex); %Useful for debuging, should be all zeros.  If there are nonzero values, increase the ammount of time the ODE solver runs for (this will slow the program down).
        %Note, longer time periods also seem to increase ODE solve step
        %intervals, which introduces its own inaccuracies so be careful
        tFin(i) = t(groundIndex);
        Ds(i) = Dp;
        rhos(i) = rhop;
        mps(i) = mp;
        times(i) = toc(runitime); %Debugging or runtime metadata
        AvgWindX(i) = mean(WindSpeedTracker(:,1)); %For non time dependent cases, WindSpeedTracker is a 1,2 vector
        AvgWindY(i) = mean(WindSpeedTracker(:,2));
        if TimeDep ~= 0
            thetas = atan2d(WindSpeedTracker(:,2),WindSpeedTracker(:,1));
            AvgTheta(i) = mean(thetas);
            Mags = sqrt(WindSpeedTracker(:,2).^2+WindSpeedTracker(:,1).^2);
            AvgMag(i) = mean(Mags);
        end
        %Print data about most recent run to command window, can eliminate for speed
        fprintf('particle %u, elapsed time %.2f, xFin = %.2f, yFin = %.2f, flight time = %.2f \n', i,times(i),xFin(i),yFin(i),tFin(i))
    end
    
    %Format Output Data to Table
    if TimeDep == 0
        OutData = table(Ds,rhos,mps,AvgWindX,AvgWindY,tFin,failedCount,times,xFin,yFin,zFin);
    else
        OutData = table(Ds,rhos,mps,AvgWindX,AvgWindY,AvgTheta,AvgMag,tFin,failedCount,times,xFin,yFin,zFin);
    end
    
    %Final Scatter Plot
    figure(5)
    scatter(xFin,yFin)
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    title('Scatterplot of Final HGM Positions')
    fprintf('%d Particles modeled with %d failures \n',numParticles,failedCount)
    %Turn off all the holds on the plots.
    %
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
    
    toc(startTime) %Elapsed time from start of script
    %patch([xPosition(1:groundIndex) nan],[yPosition(1:groundIndex) nan],[zPosition(1:groundIndex) nan],[zPosition(1:groundIndex) nan],'EdgeColor','interp','FaceColor','none')

    %function dx = createState(t,x,param)
    function dx = createState(t,x,param,TimeDep)
        %Creates the state vector containing info:
        %[Vx,Vy,Vz,ax,ay,az] from the input x = [xpos,ypos,zpos,Vx,Vy,Vz]
        Cd = param(1);
        A = param(2);
        rhof = param(3);
        mp = param(4);
        g = param(5);
        rBlower = param(6);
        modDenom = param(8);
        stdevWindChange = param(9);
        %While Particles are within blower radius, have a contant upward force.
        xPos = x(1);
        yPos = x(2);
        if xPos^2 + yPos^2 <= rBlower^2
            Fblower = param(7);
        else
            Fblower = 0;
        end

        %%%%%%%%%% THIS IS WHERE TIME DEPENDEND WIND IS CALLED %%%%%%%%%%%%%%%
        %Limits the number of times the function changes the wind speed,
        %otherwise you just get random noise.  Implemented by calling timeWind
        %every 'modDenom' seconds
        if TimeDep
            floorTime=floor(t);
            if ~mod(floorTime,modDenom) && ~ismember(floorTime,windChangeCounter) %mod(floorTime, modDenom) = 0 whever floor time is evenly divisible by modDenom
                wind = timeWind(wind,modDenom,stdevWindChange);
                windChangeCounter = [windChangeCounter,floorTime]; %Adds the current time to the windChangeCounter to ensure that the next time (which will probably be floored to the same integer) does not call timeWind again
                WindSpeedTracker = [WindSpeedTracker;wind];
            end
        end
        vwindX = wind(1);
        vwindY = wind(2);


        %Devine Velocity Components
        Vx = x(4);
        Vy = x(5);
        Vz = x(6);

        DVx = vwindX-Vx;
        DVy = vwindY-Vy;
        DVz = -Vz;
        
        %Define Acceleration Components
        %These are essentially the DEs we're solving
        ax = (Cd*A*rhof*DVx*abs(DVx)/(2*mp));
        ay = (Cd*A*rhof*DVy*abs(DVy)/(2*mp));
        az = (Cd*A*rhof*DVz*abs(DVz)/(2*mp) - g + Fblower/mp);

        %Pack derivatives into vector for ODE solver
        dx = zeros(size(x));
        dx(1) = Vx;
        dx(2) = Vy;
        dx(3) = Vz;
        dx(4) = ax;
        dx(5) = ay;
        dx(6) = az;

    end


    function  wind = timeWind(wind,modDenom,stdevWindChange)
        %time dependent wind based on a random walk
        %Change parameters to adjust random walk
        wind(1) = wind(1) + modDenom*random('Normal', 0,stdevWindChange);
        wind(2) = wind(2) + modDenom*random('Normal', 0,stdevWindChange);
        %}
    end
end