%% Manual Gamma Fitting Code

figure(1)
x = linspace(0,200e-6);
k = 4.1 %Change K and theta to fit your desired points
th = 17e-6
plot(x,gamcdf(x,k,th))
hold on
plot(x,0.1*ones(1,length(x)),'--k') %Plots lines at 10, 50, 90, 95%
plot(x,0.5*ones(1,length(x)),'--k')
plot(x,0.9*ones(1,length(x)),'--k')
plot(x,0.95*ones(1,length(x)),'--k')

perc = [0.1     0.5     0.9     0.95; 
        30e-6   65e-6   115e-6  120e-6];   %Probability Distributions entries (value,percentiles)
    
plot(perc(2,:),perc(1,:),'rx') %Plots the points you're trying to fit
str = sprintf('Fitted Gamma Cdf: k = %.2f, th = %.2e',k,th);
title(str)
hold off

%Calculates the probabilities at the specified percentile
p10=gamcdf(30e-6,k,th); %For example, this calculates what percentile 30e-6 is at (should be 0.1, but the fit isn't perfect).
p50=gamcdf(65e-6,k,th); %This number should be 50
p90=gamcdf(115e-6,k,th); %This number should be 90
p95=gamcdf(120e-6,k,th); %This number should be >=95 (hardest to fit, current fit doesns't hit)


figure(2)
x = linspace(0,200e-6);
plot(x,gampdf(x,k,th))
title('Fitted Gamma Pdf')
Ps = [p10,p50,p90,p95]*100 %Change K,th to make these fit 10, 50, 90, 95 as well as possible

sum((4/3*pi*(gamcdf(x(x>177e-6),k,th)/2).^3))/sum((4/3*pi*(gamcdf(x,k,th)/2).^3))*100 
%This calculates the percentage by mass of bubbles that are greater than 177 microns
%From Datasheet "Using a 10 gram sample on a U.S. number 80 standard sieve
%(177 microns), a maximum of five (5) percent by weight glass
%bubbles will be retained on the sieve.
%}


%% First Attempt
x0 = [k,th];
A = [0, 0];
b = [0];
param = fmincon(@minMe,x0,A,b);

%% Preliminary Gamma Code
%Plots pdfs and Cdf's of gamma across different K's and thetas to help give
%intuition about how the gamma distribution behaves.
%{
figure(1)
x = linspace(0,200e-6);
plot(x,gampdf(x,1,17e-6))
hold on
plot(x,gampdf(x,2,17e-6))
plot(x,gampdf(x,3,17e-6))
plot(x,gampdf(x,4,17e-6))
plot(x,gampdf(x,5,17e-6))
plot(x,gampdf(x,6,17e-6))
legend('1','2','3','4','5','6')
title('Shape (k) of gamma function with th = 17e-6')
hold off

figure(2)
x = linspace(0,200e-6);
plot(x,gampdf(x,4,10e-6))
hold on
plot(x,gampdf(x,4,12e-6))
plot(x,gampdf(x,4,14e-6))
plot(x,gampdf(x,4,16e-6))
plot(x,gampdf(x,4,18e-6))
plot(x,gampdf(x,4,20e-6))
legend('10','12','14','16','18','20')
title('Scale (th) of gamma function e-6 with k = 4')
hold off

figure(3)
x = linspace(0,200e-6);
plot(x,gamcdf(x,1,17e-6))
hold on
plot(x,gamcdf(x,2,17e-6))
plot(x,gamcdf(x,3,17e-6))
plot(x,gamcdf(x,4,17e-6))
plot(x,gamcdf(x,5,17e-6))
plot(x,gamcdf(x,6,17e-6))
title('Shape (k) of gamma function with th = 17e-6')
plot(x,0.1*ones(1,length(x)),'--k')
plot(x,0.5*ones(1,length(x)),'--k')
plot(x,0.9*ones(1,length(x)),'--k')
legend('1','2','3','4','5','6','10%','50%','90%')
hold off

figure(4)
x = linspace(0,200e-6);
plot(x,gamcdf(x,4,10e-6))
hold on
plot(x,gamcdf(x,4,12e-6))
plot(x,gamcdf(x,4,14e-6))
plot(x,gamcdf(x,4,16e-6))
plot(x,gamcdf(x,4,18e-6))
plot(x,gamcdf(x,4,20e-6))
title('Scale (th) of gamma function e-6 with k = 4')
plot(x,0.1*ones(1,length(x)),'--k')
plot(x,0.5*ones(1,length(x)),'--k')
plot(x,0.9*ones(1,length(x)),'--k')
legend('10','12','14','16','18','20','10%','50%','90%')
hold off

%}


%% 2 sided Gaussian Code (Archival)
%Particle Diameter Distribution
%{
Dp_mean = 65e-6;                %meters
Dp_90 = 115e-6; %90th percentile diameter
Dp_10 = 30e-6; %10th percentile diameter
Dp_max = 420e-6; %maximum particle diameter
Dp_min = 10e-6; %minimum particle diameter
Z_90 = 1.28; %Z score for 90th (and 10th) percentile
Z_10 = 1.125; %Note, not the same as tabular Z score due to min diameter cutoff
Dp_standarddevU = (Dp_90-Dp_mean)/Z_90; %Stdev above the mean
Dp_standarddevL = (Dp_mean-Dp_10)/Z_10; %Stdev below the mean
ZDmax = (Dp_max-Dp_mean)/Dp_standarddevU; %Max Z score
ZDmin = (Dp_min-Dp_mean)/Dp_standarddevL; %Min Z score
NumPart=1000000;

for i=1:NumPart
    DpOK=0;
    Dpsign = randi(2)*2-3; %randomly gives -1 or 1
    while DpOK ~=1    
        Dpstd = Dpsign*abs(random('Normal', 0, 1)); %Normally distributed standard deviation (Z score)
        if Dpsign==1 && Dpstd<ZDmax %Checks sign and whether Z is too big
            Dp=Dp_mean+Dp_standarddevU*Dpstd; %Calculates Dp from mean, stdev, and Z score
            DpOK=1; %Breaks while loop
        elseif Dpsign==-1 && Dpstd>ZDmin %Checks sign and whether Z is too small
            Dp=Dp_mean+Dp_standarddevL*Dpstd;
            DpOK=1;
        end
    end
    Dps(i)=Dp;
end

figure(1)
histogram(Dps,100)
%Checks the percentiles.  Should be 10, 10, and 5
p10 = sum(Dps<30e-6)/NumPart*100
p90 = sum(Dps>115e-6)/NumPart*100
p95 = sum(Dps>120e-6)/NumPart*100

rhop = 0.125*1000;
mass = rhop*4/3*pi*(Dps/2).^3;
mass177=sum(mass(Dps>177e-6))/sum(mass)*100;
%}


function res = minMe(x)
    P = [0.1     0.5     0.9     0.95; 
        30e-6   65e-6   115e-6  120e-6];
    
    k = x(1);
    th = x(2);
    cdfProbDiff = zeros(1,length(P));
    for i = 1:length(P)
        cdfProbDiff(i) = P(1,i) - gamcdf(P(2,i),k,th);  %difference between expected and actual percentile
    end 
    cdfProbDiff = cdfProbDiff.^2;           %square values to minimize sum.
    res = sum(cdfProbDiff);                 %sum of squared differences
end
