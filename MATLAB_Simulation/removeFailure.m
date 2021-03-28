%Strip failed runs from data

%n is the index of the values you're removing, data is the row vector
%of datasets you're changing (where each dataset is a mx1 column
%vector) Input them all as one variable

%Make sure to list outputs when running the function, or it won't work

%I'll fix the no output case and the multiple input case later if
%needed

%Example Function Call 
%fails = [1,9,28,36,75,79,85,95];
%[xFinClean,yFinClean,initWindXClean,initWindYClean,DsClean,rhosClean,timesClean]=removeFailure(fails,[xFin,yFin,initWindX,initWindY,Ds,rhos,times]);
%save('nTD_3_20.mat','xFin','yFin','initWindX','initWindY','Ds','rhos','times','xFinClean','yFinClean','initWindXClean','initWindYClean','DsClean','rhosClean','timesClean','fails')

%Particles data is printed to the command window while code is running, use
%this to identify which particles failed.
%NOTE: Error will appear before the failed particle so for case:

%particle 1 ...
%Warning: Failure
%particle 2 ...

%the failed particle is particle 2
function varargout=removeFailure(varargin)
    n = varargin{1};
    n = sort(n);
    data=varargin{2:end};

    index=0;
    for i=n
        data(i-index,:)=[];
        index=index+1; %accounts for vector size decreasing
    end
    
    for iout = 1:nargout
        varargout{iout}=data(:,iout);
    end
   
end

    
        