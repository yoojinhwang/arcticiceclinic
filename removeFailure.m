%Strip failed runs from data

%n is the index of the values you're removing, data is the row vector
%of datasets you're changing (where each dataset is a mx1 column
%vector) Input them all as one variable

%Make sure to list outputs when running the function, or it won't work

%I'll fix the no output case and the multiple input case later if
%needed

%Example Function Call 
%[xFinClean,yFinClean,timesClean]=removeFailure([67,70,71,87],[xFin,yFin,times']);
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

    
        