%Idk histograms
load('TD_4_20_PresentationData.mat')
xFin = OutData.xFin;
yFin = OutData.yFin;
figure(1)
clf
h = histogram2(xFin,yFin,[40,40],'FaceColor','flat','ShowEmptyBins','on','DisplayStyle','tile');
colorbar
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Histogram of Final HGM Positions')

xplotvec = linspace(0,3400);
figure(2)
clf
scatter(xFin,yFin)
hold on
plot(xplotvec,500/2000.*xplotvec,'--r')
plot(xplotvec,-500/2000.*xplotvec,'--r')

xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Scatterplot of Final HGM Positions')