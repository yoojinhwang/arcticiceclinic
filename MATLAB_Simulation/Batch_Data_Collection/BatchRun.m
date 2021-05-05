%Batch Run
numParticles = 300; %Number of particles to run, integer
TimeDep = 1; % 0 for non time dependent, otherwise, wind varies with time

%Time Dependent Variables (only used if TimeDep ~= 0)
WindChangeFreq = 2; %frequency of wind change (s), integer
stdevWindChange = 0.5; %Standard Deviation of wind change, float

% Change the filename to the desired .mat file
filename = 'TD_4_17_HighCDBen2.mat';% 'TD_4_14.mat'
%Keeps you from overwriting files
if isfile(filename)
    disp('File already exists, if you wish to overwrite, please delete the file')
    return
elseif length(strsplit(filename,'.'))==1
    noExt = input('Filename has no extension, by default, a .mat file will be created.\nPress return/enter without any other input to cancel.\nPress any other key and enter to continue:');
    if isempty(noExt)
        return
    end
end


OutData = ParticleModel(numParticles,TimeDep,WindChangeFreq,stdevWindChange);

save(filename,'OutData')


