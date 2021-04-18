%stats
%load data or use data preloaded
%{
load('TD_3_27_1.mat')
mass = rhos.*4/3*pi.*(Ds/2).^3;
TMass = table(rhos,Ds,mass,xFin);
TNoMass = table(rhos,Ds,xFin);
mdlMass = fitlm(TMass)
mdlNoMass = fitlm(TNoMass)
%}
load('TD_4_17_HighCDBen.mat')
mass = OutData.mps;
Ds = OutData.Ds;
rhos = OutData.rhos;
xFin = OutData.xFin;
AvgWindX = OutData.AvgWindX;
yFin = OutData.yFin;
mdlTNoMass = table(Ds,rhos,AvgWindX,xFin);
mdlTMass = table(Ds,rhos,mass,AvgWindX,xFin);
mdlTJustMass = table(mass,AvgWindX,xFin);
mdlTNoWind = table(Ds,rhos,mass,xFin);
mdlNoMass = fitlm(mdlTNoMass)
mdlMass = fitlm(mdlTMass)
mdlJustMass = fitlm(mdlTJustMass)
mdlNoWind = fitlm(mdlTNoWind)
AvgX = mean(xFin)
StDevX = std(xFin)
AvgY = mean(yFin)
StDevY = std(yFin)