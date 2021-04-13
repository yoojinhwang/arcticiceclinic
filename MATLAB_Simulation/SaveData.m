%Save data to file
load nTD_3_27_1.mat
filename = 'IceData_3_27_1.xlsx';
sh = 'nTD';

%fails = fails';
Dist = sqrt(xFin.^2+yFin.^2);
mass = rhos.*4/3*pi.*(Ds/2).^3;
T = table(rhos,Ds,mass,initWindX,initWindY,times,xFin,yFin,Dist);

writetable(T,filename,'Sheet',sh)






%Tclean = table(xFinClean,yFinClean,rhosClean,DsClean,initWindXClean,initWindYClean,timesClean);
%Tfails = table(fails);
%writetable(Tclean,filename,'Sheet',sh,'Range','H1')
%writetable(Tfails,filename,'Sheet',sh,'Range','O1')