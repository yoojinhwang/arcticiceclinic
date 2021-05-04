function OutData = ParticleModel(numParticles,solType,Twind,stdevWindChange,graph)
    %numParticles:  Number of iterations/particles
    %solType:   The solution type for the model:
    %           "ANA" -- Anaytical Solution"
    %           "CW" -- Constant Wind Solution"
    %           "TDW" -- Time Dependent Wind Solution"
    %           Must be String enclosed by " " 
    %Twind:  Only relevant in TD models, frequency of wind change in
    %        seconds.  Must be an integer
    %stdevWindChange:   Only relevant in TD models, standard deviation of
    %                   random wind speed changes (gaussian distribution)
    %graph: Boolean enabling graphs, set to zero to significantly improve
    %       runtime
    
    %
    %OutData:   A table of output data.  Can be configured to fit your data
    %           needs

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

    %Other
    Cd = 0.659692; %FURTHER WORK: Vary CD with reynolds number
    g = 9.81;

    %Adjustment to ode solver to prevent squiggles, increases comp time
    options = odeset('RelTol',1e-5,'AbsTol',1e-7); 
    
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
    times = zeros(numParticles,1);
    failed= zeros(numParticles,1);

    for i = 1:numParticles
        %%Create Random Parameters
        runitime=tic;
        %Random Wind
        vwind = random('Normal', Vwind_mean,Vwind_standarddev);
        thetaWind = random('Normal',thetaWind_mean,thetaWind_standarddev);
        wind(1) = vwind*cos(thetaWind); %vwindX
        wind(2) = vwind*sin(thetaWind); %vwindY
        WindSpeedTracker = wind; %resets WindSpeedTracker variable and loads first wind speed


        %Random Diameter
        DpOK=0;
        while DpOK ~=1    
            Dp =random('Gamma', k, th);
            %Dp = random('Normal',Dp_mean,Dp_standarddev); %Use for standard gaussian distribution
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
        
        if solType == "ANA"
            %Analytical Solution
            %Since the code only solves for final position, plots are
            %irrelevant
            graph = 0;

            %Initial Conditions
            initX = 0;
            initY = 0;
            initZ = htanker;
            initVx = 0;
            initVy = 0;
            initVz = vblower;
            %Define wind components
            vwindX = wind(1);
            vwindY = wind(2);

            %Coeffecients for the z equation
            a = mp;
            b = (Cd*A*rhof)/(2);
            c = mp*g;
            c1 = 1/sqrt(b*c)*atanh(-initVz*sqrt(b/c));
            c2 = initZ + a/b*log(cosh(c1*sqrt(b*c)));

            %Actual Equation:t_z0 = a/sqrt(b*c)*acosh(exp(b/a*c2)) - a*c1;
            %Approximate: acosh(exp(x)) = x + 0.693
            t_z0 = a/sqrt(b*c)*(b/a*c2 + 0.693) - a*c1;
            %Accounting for initVz condition being outside the real domain of atanh
            t_z0 = real(t_z0);  
            tFin(i) = t_z0;

            %Coefficient for x, and y equations
            a = Cd*rhof*A/(2*mp);           
            %x-equation for the analytical solution
            c1x = 1/(vwindX - initVx);
            c2x = initX - 1/a*(vwindX*c1x - log(c1x));
            xFin(i) = 1/a*(vwindX*(a*t_z0 + c1x) - log(a*t_z0 + c1x)) + c2x;

            %The form of the y-equation depends on if vwindY is positive or negative
            if vwindY >= 0             
                c1y =  1/(vwindY - initVy);
                c2y = initY - 1/a*(vwindY*c1y - log(c1y));
                yFin(i) = 1/a*(vwindY*(a*t_z0 + c1y) - log(a*t_z0 + c1y)) + c2y;
            else
                c1y = 1/(vwindY*sqrt(2))*atanh((vwindY - initVy)/(vwindY*sqrt(2)));
                c2y = initY + 1/a*log(cosh(vwindY*sqrt(2)*c1y));
                yFin(i) = c2y + vwindY*t_z0 - 1/a*((vwindY*sqrt(2)*c1y) - 0.693);
                %Actual Equation:c2y + vwindY*t_z0 - 1/a*log(cosh(vwindY*sqrt(2)*c1y)
                %Appriximate: acosh(exp(x)) = x + 0.693
            end
            %End Analytical Solution
        elseif solType == "TDW" || solType == "CW"
            if solType == "TDW"
                TimeDep = 1;
            else
                TimeDep = 0;
            end
            %Pack parameters into vector for createState Function
            parameters = zeros(6,1);        %Adjust if you want to add more parameters
            parameters(1) = Cd;
            parameters(2) = A;
            parameters(3) = rhof;
            parameters(4) = mp;
            parameters(5) = g;
            parameters(6) = rBlower;
            parameters(7) = FBlower;
            parameters(8) = Twind;
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
            %This may change when the mass distribution is corrected
            if Dp>50e-6
                ODETimeLim = 120; %120 and 200 seconds is a good start for K1 HGMs
            else
                ODETimeLim = 200;
            end

            %Call Ode23 to solve governing equations
            windChangeCounter = 0; %Ensures wind changes just once per time interval
            [t,y] = ode45(@(t,x)createState(t,x,parameters,TimeDep),[0,ODETimeLim],x0,options);
            xPosition = y(:,1);
            yPosition = y(:,2);
            zPosition = y(:,3);
            xVelocity = y(:,4);
            yVelocity = y(:,5);
            zVelocity = y(:,6);

            %Find the index for when particles reach ground.
            [minZ,groundIndex] = min(abs(zPosition));
            %Check to see if they actually hit the ground or just ran out of time

            if minZ>0.5 && t(groundIndex)==ODETimeLim %This is an estimation, not an actual count
                fprintf('Error: Particle may not have hit ground, Zpos = %.2f, t = %f\n',minZ,t(groundIndex))
                failed(i) = 1;
            elseif minZ>0.5 %Adjust this value to modify final height tolerance
                fprintf('Warning: Final particle height outside of tolerance\n')
            end

        else
            typeSol = [ "'ANA' -- Anaytical Solution";
                        "'CW' -- Constant Wind Solution";
                        "'TDW' -- Time Dependent Wind Solution"];
            fprintf('Solution Type not recognized. Please choose one of the following for SolType:\n')
            graph = 0;
            disp(typeSol)
        end
        if graph
            figure(1)
            subplot(2,1,1)
            plot(t(1:groundIndex), zVelocity(1:groundIndex));
            hold on
            title('Z Velocity and Position')
            xlabel('time (s)')
            ylabel('Z Velocity (m/s)')
            subplot(2,1,2)
            plot(t(1:groundIndex), zPosition(1:groundIndex));
            hold on
            xlabel('time (s)')
            ylabel('Z Position (m/s)')

            figure(2)
            subplot(2,1,1)
            plot(t(1:groundIndex), xVelocity(1:groundIndex));
            hold on
            title('X Velocity and Position')
            xlabel('time (s)')
            ylabel('X Velocity (m/s)')
            subplot(2,1,2)
            plot(t(1:groundIndex), xPosition(1:groundIndex));
            hold on
            xlabel('time (s)')
            ylabel('X Position (m/s)')

            figure(3)
            subplot(2,1,1)
            plot(t(1:groundIndex), yVelocity(1:groundIndex));
            hold on
            title('Y Velocity and Position')
            xlabel('time (s)')
            ylabel('Y Velocity (m/s)')
            subplot(2,1,2)
            plot(t(1:groundIndex), yPosition(1:groundIndex));
            hold on
            xlabel('time (s)')
            ylabel('Y Position (m/s)')

            figure(4)
            plot(xPosition(1:groundIndex), yPosition(1:groundIndex))
            hold on
            xlabel('X Position')
            ylabel('Y Position')
            title('Particle Pathlines')

            figure(5)
            plot3(xPosition(1:groundIndex),yPosition(1:groundIndex),zPosition(1:groundIndex),'LineWidth',3) 
            title('Particle Pathlines')
            xlabel('X Position')
            ylabel('Y Position')
            zlabel('Z Position')
            grid on
            hold on
        end
        
        
        %Output Data Collection
        if solType ~= "ANA"
            xFin(i) = xPosition(groundIndex);
            yFin(i) = yPosition(groundIndex);
            zFin(i) = zPosition(groundIndex);
            tFin(i) = t(groundIndex);
        end
        Ds(i) = Dp;
        rhos(i) = rhop;
        mps(i) = mp;
        times(i) = toc(runitime); %Debugging or runtime metadata
        AvgWindX(i) = mean(WindSpeedTracker(:,1)); %For non time dependent cases, WindSpeedTracker is a 1,2 vector
        AvgWindY(i) = mean(WindSpeedTracker(:,2));
        if solType == "TDW"
            thetas = atan2d(WindSpeedTracker(:,2),WindSpeedTracker(:,1));
            AvgTheta(i) = mean(thetas);
            Mags = sqrt(WindSpeedTracker(:,2).^2+WindSpeedTracker(:,1).^2);
            AvgMag(i) = mean(Mags);
        end
        %Print data about most recent run to command window, can comment for speed or clarity
        fprintf('particle %u, elapsed time %.2f, xFin = %.2f, yFin = %.2f, flight time = %.2f \n', i,times(i),xFin(i),yFin(i),tFin(i))
    end
    
    %Format Output Data to Table
    if solType == "CW"
        OutData = table(Ds,rhos,mps,AvgWindX,AvgWindY,tFin,failed,times,xFin,yFin,zFin);
    elseif solType == "TDW"
        OutData = table(Ds,rhos,mps,AvgWindX,AvgWindY,AvgTheta,AvgMag,tFin,failed,times,xFin,yFin,zFin);
    elseif solType == "ANA"
        OutData = table(Ds,rhos,mps,AvgWindX,AvgWindY,times,tFin,xFin,yFin);
    else
        typeSol = [ "'ANA' -- Anaytical Solution";
                    "'CW' -- Constant Wind Solution";
                    "'TDW' -- Time Dependent Wind Solution"];
        fprintf('Solution Type not recognized. Please choose one of the following for SolType:\n')
        disp(typeSol)
    end
    
    fprintf('%d Particles modeled with %d failures \n',numParticles,sum(failed))
    
    %Final Scatter Plot
        figure(6)
        scatter(xFin,yFin)
        xlabel('X Position (m)')
        ylabel('Y Position (m)')
        title('Scatterplot of Final HGM Positions')

    if graph
        %Turn off all the holds on the plots.
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

        figure(5)
        hold off
    end
    %}
    
    toc(startTime) %Elapsed time from start of script
    
    function dx = createState(t,x,param,TimeDep)
        %Creates the state vector containing info:
        %[Vx,Vy,Vz,ax,ay,az] from the input x = [xpos,ypos,zpos,Vx,Vy,Vz]
        Cd = param(1);
        A = param(2);
        rhof = param(3);
        mp = param(4);
        g = param(5);
        rBlower = param(6);
        Twind = param(8);
        stdevWindChange = param(9);
        xPos = x(1);
        yPos = x(2);
        %While Particles are within blower radius, have a contant upward force.
        if xPos^2 + yPos^2 <= rBlower^2
            Fblower = param(7);
        else
            Fblower = 0;
        end

        %Call Time Dependent Wind
        if TimeDep
            floorTime=floor(t);
            if ~mod(floorTime,Twind) && ~ismember(floorTime,windChangeCounter) %mod(floorTime, Twind) = 0 when floor time is evenly divisible by Twind
                wind = timeWind(wind,Twind,stdevWindChange);
                windChangeCounter = [windChangeCounter,floorTime]; %Adds the current time to the windChangeCounter to ensure that the next time step (which will probably be floored to the same integer) does not call timeWind again
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


    function  wind = timeWind(wind,Twind,stdevWindChange)
        %time dependent wind based on a random walk
        wind(1) = wind(1) + Twind*random('Normal', 0,stdevWindChange);
        wind(2) = wind(2) + Twind*random('Normal', 0,stdevWindChange);
    end
end