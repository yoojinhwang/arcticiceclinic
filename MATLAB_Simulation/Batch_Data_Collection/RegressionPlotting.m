%Plotting Regression Results
load('NTD_4_20_PresentationData.mat')
mps = OutData.mps;
xFin = OutData.xFin;
yFin = OutData.yFin;
Ds = OutData.Ds;
rhos = OutData.rhos;
AvgWindX = OutData.AvgWindX;
mdlTJustMass = table(mass,AvgWindX,xFin);
mdlJustMass = fitlm(mdlTJustMass)
bJustMass = mdlJustMass.Coefficients.Estimate;

%{
figure(1)
scatter3(Ds,rhos,xFin,'Filled')
title('Regression on Ds and rho')
xlabel('Ds')
ylabel('rhos')
zlabel('xFin')
hold on
DsFit = min(Ds):10^-5:max(Ds);
rhosFit = min(rhos):4:max(rhos);
[DSFIT,RHOSFIT]=meshgrid(DsFit,rhosFit);
YFIT = bNoMass(1)+bNoMass(2)*DSFIT+bNoMass(3)*RHOSFIT; %Regression Coefs, if you don't have them, run them yourself
mesh(DSFIT,RHOSFIT,YFIT)
hold off


figure(2)
scatter3(Ds,rhos,xFin,'Filled')
title('Regression on Ds, rho, and Mass')
xlabel('Ds')
ylabel('rhos')
zlabel('xFin')
hold on
DsFit = min(Ds):10^-5:max(Ds);
rhosFit = min(rhos):4:max(rhos);
[DSFIT,RHOSFIT]=meshgrid(DsFit,rhosFit);
YFIT = bMass(1)+bMass(2)*DSFIT+bMass(3)*RHOSFIT+bMass(4)*(4*pi/24)*RHOSFIT.*DSFIT.^3; %Regression Coefs, if you don't have them, run them yourself
mesh(DSFIT,RHOSFIT,YFIT)
hold off
%}

figure(3)
scatter3(mps,AvgWindX,xFin,'Filled')
title('Regression on Mass and Wind Velocity')
xlabel('HGM Mass')
ylabel('Average X Wind Velocity')
zlabel('Final X Position')
hold on
MassFit = linspace(min(mps),max(mps),10);
VelFit = linspace(min(AvgWindX),max(AvgWindX),10);
[MASSFIT,VELFIT]=meshgrid(MassFit,VelFit);
YFIT = bJustMass(1)+bJustMass(2)*MASSFIT+bJustMass(3)*VELFIT; %Regression Coefs, if you don't have them, run them yourself
mesh(MASSFIT,VELFIT,YFIT,'FaceAlpha','0')
hold off

mdlTMass = table(Ds,rhos,mass,AvgWindX,xFin);
mdlMass = fitlm(mdlTMass)
AvgX = mean(xFin);
StDevX = std(xFin);
AvgY = mean(yFin);
StDevY = std(yFin);
Avgs = [AvgX;AvgY];
Stds = [StDevX;StDevY];
NTDTab = table(Avgs,Stds,'RowNames',{'X Position','Y Position'})