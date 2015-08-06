% estSTA_Data_3_1 This "function" evalues the real linear Kernel and checks the
% goodness of fit to a predefined model 

% Load data
D = LoadDataIntoProcess( 1 ); % loads Segev's data
process = D.Data;

% define upper and lower bounderies for optimization
LowConstr = [eps -100 -100 0]; % [sig f mu phi]
UppConstr = [10 100 100 pi]; % [sig f mu phi]

dt_STIM = process.stimulus.dt;
Fs_STIM = 1/dt_STIM;
time_STIM = linspace(0,  length(process.stimulus.Yuncorr)*dt_STIM,  size(process.stimulus.Yuncorr,1));
STA_LENGTH = 35;
process = process.CalcSTA( STA_LENGTH, 5);

warning('off');

%         figure(1);
%         plot(time,targetfun(time));
%         figure(2);plot(time,modelfun(p(1),p(2),p(3),p(4)))
%         figure(1);plot(time,Y_est)
%         pause(1);
%% Optimization 

% using the same optimization procedure as in section 2.1
x0 = [1 1 1 1]; % start point
time = process.STA.timeSTA;
numcoeff = 4;
N = length(time);
dfe = N - numcoeff; % degrees of freedom - number of observations - num of coefficients
realSTA = {process.STA.realSTA}';
% modelfun = @(sig, f, mu, phi) exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ); % original Kernel model
% costfun = @(Y) (@(x) sum( (Y' - modelfun(x(1),x(2),x(3),x(4))).^2 )); % SA with direct search  [sig f mu phi]
modelfun = @(sig, f, mu, phi) normax( exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ) ); % with normax
costfun = @(Y) (@(x) sum( ( Y' - modelfun(x(1),x(2),x(3),x(4))).^2 )); % SA with direct search  [sig f mu phi]

% [p_est, output, goodness, res, est] = deal(cell(process.numChannels,1));
estSTA = cell(process.numChannels,1);
%% Solving with a parfor Loop

[estSTA, res] = deal(zeros(process.numChannels, STA_LENGTH));
% goodness = mat2cell(zeros(process.numChannels, 1));
goodness = cellfun(@(x) struct('sse',[],'rsquare',0,'dfe',[],'adjrsquare',[],'rmse',[]), num2cell(1:process.numChannels));
processArray = (1:process.numChannels); % the processes which needs to be estimated

% parameters for mulitple runs of the optimization
AreAllOptim = true( process.numChannels, 1); % flag if the fitting has R2 > .85
maxRun = 20; % maximum nubmbers of runs of the optimization
numIterGlobal = 0; % index for iteration number of whole process

optimData = struct('LowConstr',LowConstr,'UppConstr',UppConstr,'numcoeff',...
    numcoeff,'process',process,'STA_LENGTH',STA_LENGTH,'costfun',costfun,...
    'modelfun',modelfun,'N',N,'dfe',dfe);
[ estSTA, res, goodness, estPARAM ] = multiOptim_PARFOR( realSTA, optimData );

range = (1:process.numChannels);
optimResults = struct('time',time,'data',realSTA,'estimationSTA',estSTA,...
    'estimationPARAM',estPARAM,'Stat_goondness',goodness);               
%% Parralel loop through all STAs
% 2DO - each worker overrides 2 calculations (except of the last one)
% generalize to the case where the local parrallel profile has less workers
% than data

processArray = (1:process.numChannels); % the processes which needs to be estimated
numIterD = distributed( processArray );
% numIterD = distributed( (1:6) );

% parameters for mulitple runs of the optimization
AreAllOptim = true( process.numChannels, 1); % flag if the fitting has R2 > .85
maxRun = 3; % maximum nubmbers of runsof the optimization
numIterGlobal = 0; % index for iteration number of whole process
optimResults = cell(maxRun,1);

while any( AreAllOptim ) && numIterGlobal < maxRun % if R^2 > .85 and didn't exceed maxRun
    spmd
    %         [estSTA, res] = deal(codistributed.zeros(process.numChannels,STA_LENGTH));
    %         goodness =  codistributed.cell(process.numChannels,1);
            numIterLP = getLocalPart ( numIterD );
            for n = numIterLP
    %         for n = numIterLP(1):numIterLP(end)      
                fprintf('working on %d''th STA fitting ...\n', n);    

                costIter = costfun(normax( realSTA{n} ));
    %             [p_est,~,~,output] = annealData2Acorr(x0, LowConstr,UppConstr,costIter);          
    %             [p_est,~,~,output] = gps_Data2STA(x0 , LowConstr,UppConstr, costIter);
                [p_est,~,~,output] = ga_Data2STA( numcoeff, LowConstr,UppConstr, costIter);
    %             fprintf('current n is %d \n' ,n); 
                estSTA = modelfun( p_est(1), p_est(2), p_est(3), p_est(4)); % [sig f mu phi]   
    %             plotFit( time, normax(realSTA{n}) , estSTA, 'STA 1 fit')

                % Evaluating goodness of optimization and fit
                res = normax(realSTA{n}) - estSTA';  
                goodness = iGoodnessStructure(estSTA,[],res,dfe,N);
                
                if goodness.rsquare < .85
                    AreAllOptim(n) = false;
                end
    %             estSTA(n,:) = modelfun( p_est(1), p_est(2), p_est(3), p_est(4)); % [sig f mu phi]   
    %             res(n,:) = normax(realSTA{n}) - estSTA';  
    %             goodness(n) = iGoodnessStructure(estSTA,[],res,dfe,N);
            end
    end  
    AreAllOptim = logical( prod(cell2mat(AreAllOptim(1:end)),2) ); % uniting the composite values to 
    disp(['Are All Optimized? ',num2str( AreAllOptim' )])
    numIterGlobal = numIterGlobal + 1;
    numIterD = distributed(  processArray(AreAllOptim) );
    optimResults{numIterGlobal} = struct('time',time,'data',{realSTA(range)},'estimation', {estSTA(range)},'Stat_goondness',{goodness(range)},'optimOut', {output(range)});               
end
range = (1:process.numChannels);
optimResults = struct('time',time,'data',{realSTA(range)},'estimation', {estSTA(range)},'Stat_goondness',{goodness(range)},'optimOut', {output(range)});               

temp = goodness(range);cellfun(@(x)disp(x.rsquare),temp)
