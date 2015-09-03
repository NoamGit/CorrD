function [ lambda_corr ] = ppcorr(obj, bias, correctR0 , method)

% R=ppcorr(obj.spiketimes, obj.dt, obj.obj.maxlags, bias, correctR0) calculates the correlation
%       structure of the point processes
% obj.spiketimes - cell array (may be a vector ) of spike obj.spiketimes of different units vector (for a single process) or
% obj.dt - desired temporal resolution
% maxlags - maxlags (as in xcorr) - full length if not specified
% bias - bias flag (as in xcorr, on of: {'none', 'biased', ['unbiased']})
% correctR0 - flag whether to remove a 'delta' at Rxx(tau=0)
% (1 - correct (MGF), 2 - remove, 0 - don't remove). Default set to 0.

% check the input apply defaults
if nargin<2
    obj.maxlags=[];
end

if nargin<3
    bias='unbiased';
end

if nargin<4
    correctR0=0;
end

% calculate the correlation structure
nCells=length(obj.spiketimes);
lambda_corr=cell(nCells);

% a patch for a single vector obj.spiketimes input
switch method  
    case 'cross'
        vFlag = false;
    otherwise 
        vFlag = true;
        lambda_Acorr=cell(nCells,1);
end

for iCell=1:nCells
    obj.spiketimes{iCell}=obj.spiketimes{iCell}-min(obj.spiketimes{iCell});
    Tmax=max(obj.spiketimes{iCell});
    Tmin=min(obj.spiketimes{iCell});
    counts1=histc(obj.spiketimes{iCell}, Tmin:obj.dt:Tmax);
    if isempty(obj.maxlags)
        obj.maxlags=length(counts1);
    end
    lambda_corr{iCell, iCell}=xcorr(counts1, obj.maxlags, bias)/obj.dt^2;

    % correct delta at tau==0
    if correctR0 == 1
        % this implementation is using the moment generating function
        zerolag=(length(lambda_corr{iCell, iCell})+1)/2;
        E=ppmean(obj.spiketimes{iCell});
        lambda_corr{iCell, iCell}(zerolag)=lambda_corr{iCell, iCell}(zerolag)-E/obj.dt; 
    elseif correctR0 == 2
        zerolag=(length(lambda_corr{iCell, iCell})+1)/2;
        lambda_corr{iCell, iCell}(zerolag)=[]; 
    end
    
    if vFlag % calculates only AC (not Cross-Corr) if flag is up
        lambda_Acorr{iCell} = lambda_corr{iCell, iCell};
        continue;
    end
    
    for jCell=iCell+1:nCells
        Tmax=max(max(obj.spiketimes{iCell}), max(obj.spiketimes{jCell}));
        Tmin=min(min(obj.spiketimes{iCell}), min(obj.spiketimes{jCell}));
        counts1=histc(obj.spiketimes{iCell}, Tmin:obj.dt:Tmax);
        counts2=histc(obj.spiketimes{jCell}, Tmin:obj.dt:Tmax);
        if isempty(obj.maxlags)
            obj.maxlags=length(counts1);
        end
        lambda_corr{iCell, jCell}=xcorr(counts1, counts2, obj.maxlags, bias)/obj.dt^2;
        % using symmetry to save computations
        lambda_corr{jCell, iCell}=flipud(lambda_corr{iCell, jCell});
    end
end

% return a vector (and not a cell array) if obj.spiketimes was a vector
if vFlag
    lambda_corr = lambda_Acorr;
end

