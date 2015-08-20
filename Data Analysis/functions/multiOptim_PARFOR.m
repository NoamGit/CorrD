function [ estSTA,res, goodness, estPARAM ] = multiOptim_PARFOR( data, paramSTRUCT )
%[ estSTA,res, goodness ] = multiOptim_PARFOR( paramSTRUCT ) STA fitting
% subfunction. defined for order and debugging option

LowConstr = paramSTRUCT.LowConstr;
UppConstr = paramSTRUCT.UppConstr;
numcoeff = paramSTRUCT.numcoeff;
process = paramSTRUCT.process;
STA_LENGTH = paramSTRUCT.STA_LENGTH;
costfun = paramSTRUCT.costfun;
modelfun = paramSTRUCT.modelfun;
N =  paramSTRUCT.N;
dfe =  paramSTRUCT.dfe;

[estSTA, res] = deal(zeros(process.numChannels, STA_LENGTH));
estPARAM = zeors(process.numChannels, numcoeff);
goodness = cellfun(@(x) struct('sse',[],'rsquare',0,'dfe',[],'adjrsquare',[],'rmse',[]), num2cell(1:process.numChannels));
processArray = (1:process.numChannels); % the processes which needs to be estimated

% parameters for mulitple runs of the optimization
maxRun = 20; % maximum nubmbers of runsof the optimization
numIterGlobal = 0; % index for iteration number of whole process

while numIterGlobal < maxRun % if R^2 > .85 and didn't exceed maxRun  
    parfor indx = processArray
        costIter = costfun(normax( data{indx} ));
        p_est_ITER = ga_Data2STA( numcoeff, LowConstr,UppConstr,costIter );
        
        estSTA_ITER = modelfun( p_est_ITER(1), p_est_ITER(2), p_est_ITER(3), p_est_ITER(4) )'; % [sig f mu phi] 
        res_ITER = normax(data{indx}) - estSTA_ITER;  
        goodness_ITER = iGoodnessStructure(estSTA_ITER,[],res_ITER ,dfe,N);
        if goodness_ITER.rsquare > goodness(indx).rsquare
            estSTA(indx,:) = estSTA_ITER;
            res(indx,:) = res_ITER;
            goodness(indx) = goodness_ITER;
            estPARAM(indx) = p_est_ITER;
        end
    end
    numIterGlobal = numIterGlobal + 1;
end


