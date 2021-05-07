%Plotting Regression Results
load('TD_4_20_PresentationData.mat')
mps = OutData.mps;
xFin = OutData.xFin;
Ds = OutData.Ds;
D2 = Ds.^2;
rhos = OutData.rhos;
AvgWindX = OutData.AvgWindX;
mdlT = table(Ds,D2,rhos,AvgWindX,xFin);
mdl = fitlm(mdlT)
b = mdl.Coefficients.Estimate
avgRho = 0.125*1000;

figure(1)
scatter3(Ds,AvgWindX,xFin,'Filled')
title('Time Dependent Wind Regression')
xlabel('Ds')
ylabel('Average X Wind')
zlabel('Final X Position')
hold on
DsFit = linspace(min(Ds),max(Ds),10);
rhosFit = linspace(min(AvgWindX),max(AvgWindX),10);
[DSFIT,WINDFIT]=meshgrid(DsFit,rhosFit);
YFIT = b(1) + b(2)*DSFIT + b(5)*WINDFIT + b(3)*DSFIT.^2+avgRho*b(4);
mesh(DSFIT,WINDFIT,YFIT,'FaceAlpha','0.3')
hold off

%{
mdlTMass = table(Ds,rhos,mass,AvgWindX,xFin);
mdlMass = fitlm(mdlTMass)
AvgX = mean(xFin);
StDevX = std(xFin);
AvgY = mean(yFin);
StDevY = std(yFin);
Avgs = [AvgX;AvgY];
Stds = [StDevX;StDevY];
NTDTab = table(Avgs,Stds,'RowNames',{'X Position','Y Position'})
%}