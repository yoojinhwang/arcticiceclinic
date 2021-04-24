%Batch run to scan across frequencies/intensities
numParticles = 300; %Number of particles to run, integer
TimeDep = 1; % 0 for non time dependent, otherwise, wind varies with time

%Time Dependent Variables (only used if TimeDep ~= 0)
Freqs = [2:2:10,30]; %frequency of wind change (s), integer
Mags = [0.05,linspace(0.1,1,10),2]; %Standard Deviation of wind change, float

filename = 'TD_Parameter_Sweep_4_24';% 'TD_4_14.mat'

MeanXs = zeros(length(Freqs),length(Mags));
StdXs  = zeros(length(Freqs),length(Mags));
MeanYs = zeros(length(Freqs),length(Mags));
StdYs  = zeros(length(Freqs),length(Mags));
fails = zeros(length(Freqs),length(Mags));

for i = 1:length(Freqs)
    for j = 1:length(Mags)
        OutData = ParticleModelBen(numParticles,TimeDep,Freqs(i),Mags(j));
        MeanXs(i,j) = mean(OutData.xFin);
        StdXs(i,j) = std(OutData.xFin);
        MeanYs(i,j) = mean(OutData.yFin); 
        StdYs(i,j) = std(OutData.yFin);
        fails(i,j) = OutData.failedCount;
    end
end

save(filename,'MeanXs','StdXs','MeanYs','StdYs','Mags', 'Freqs')
        