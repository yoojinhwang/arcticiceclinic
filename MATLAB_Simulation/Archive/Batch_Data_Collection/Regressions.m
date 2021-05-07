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
load('NTD_4_20_PresentationData.mat')
mass = OutData.mps;
Ds = OutData.Ds;
D2 = Ds.^2;
D3 = Ds.^3;
rhos = OutData.rhos;
xFin = OutData.xFin;
AvgWindX = OutData.AvgWindX;
yFin = OutData.yFin;
mdlTNoMass = table(Ds,rhos,AvgWindX,xFin);
mdlTMass = table(Ds,rhos,mass,AvgWindX,xFin);
mdlTJustMass = table(mass,AvgWindX,xFin);
mdlTNoWind = table(Ds,rhos,mass,xFin);
mdlTNoMassD2 = table(Ds,D2,rhos,AvgWindX,xFin);
mdlTMassD2 = table(Ds,D2,rhos,mass,AvgWindX,xFin);
mdlTNoWindD2 = table(Ds,D2,rhos,mass,xFin);
mdlTJustD2 = table(D2,rhos,AvgWindX,xFin);
mdlTJustD2mass = table(D2,rhos,mass,AvgWindX,xFin);
mdlTNoMassD3 = table(Ds,D2,D3,rhos,AvgWindX,xFin);

mdlNoMass = fitlm(mdlTNoMass);
mdlMass = fitlm(mdlTMass);
mdlJustMass = fitlm(mdlTJustMass);
mdlNoWind = fitlm(mdlTNoWind);
mdlNoMassD2 = fitlm(mdlTNoMassD2)
mdlJustD2 = fitlm(mdlTJustD2)
mdlNoMassD3 = fitlm(mdlTNoMassD3);


%Regressions with both D2 and Mass return NaN for mass (D2 too closely
%correlated, matrix isn't full rank)

%mdlMass.Rsquared
%mdlNoMassD2.Rsquared
%mdlMassD2.Rsquared


AvgX = mean(xFin)
StDevX = std(xFin)
AvgY = mean(yFin)
StDevY = std(yFin)