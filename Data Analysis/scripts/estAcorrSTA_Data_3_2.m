% estAcorrSTA_Data_3_2 

%% load data module
% comparing parametric model to fit data in forward and backward
% procedure. The comparison is done in the nl and the CGP plane

% Loads data (also calculates CP)
    nlType = 'poly';
    D = LoadDataIntoProcess( 'Segev''s Data', nlType ); % loads Segev's data
    process = D.Data;

% Define structure of assumed linear kernel
    % we take this width because we want max lags have 0 reminder with both sampling rates
    maxlags = ceil((0.621 * sqrt(2*log(10)) * 1)/process.dt); 
    process.maxlags = maxlags;
    
%% pre-processing data module
% preproccessing and calculating properties of model : CGP & lambda acorr,
% CGP, STA ,nl

    poly_order = 2; % order of polynum to be fit to the nl    
    showFlag = 0; % compare RdN flag
 
    % DoubletCancel - change CP behaviour to be more poissionian like
    flagDoublet = 1;
    if flagDoublet
        obj = DoubletCancel(process);
    end
    
    % CalcSTA - finds the sta and potentially upsamples the stimulus
    process = process.CalcSTA('original stimulus'); 
    
    % nlestimation - estimates nl and calculates CGP 
    amp = 1;
    process = process.nlestimation( poly_order, showFlag, amp, nlType,'original stimulus', 'resample CP' ); 
        % process = process.CalcSTA('resampled stimulus'); % finds the sta and upsamples the stimulus
        % process = process.nlestimation(poly_order, showFlag, 'resample stimulus'); 
        
    % calculates acorr and can calc CGP ( showFlag 1 yields good samples 3 6 10 11 13 14 17 18 21 22 23)
    flagCompare = 0; % Compare = 1
    process = process.PreProcessCP( nlType, flagDoublet, flagCompare );
%% data visualization module
warning('off');

% Show STA estimation - results : exclude cell 16  
    [STA_trim,time_trim] = deal(cell(process.numChannels,1));
    for k = 1:process.numChannels
        idx = process.STA(k).timeSTA <= 0;
        STA_trim{k} = process.STA(k).STA(idx);
        time_trim{k} = process.STA(k).timeSTA(idx);
    end
    plotProperties = struct('time',time_trim,'method',...
        0,'num2disp',6,'title','STA cell - ','xlabel','time [sec]');
    plotCompare( STA_trim, plotProperties );

% show nl estimation - results 'exp' and 'poly' : exclude cells 1 7 9 10 13 14 16 21 22 
    showFlag = 1;
    process = process.nlestimation( poly_order, showFlag, amp, nlType,'original stimulus' ); 

% ** optional - see function description - result : exlude cells 9 
% current problem - COMPARE SCALE WITH ONENOTE
    compareInNlPlane( process );

% ** optional - see func description
% result - good statistical prop in cell 1 2 4 5 7 8 9 14 15 20 21 22  
    showISI( process );

% ** optional - 
    burstDetection( process );
%% back to data cleansing


%% Data fitting and optimization module

% define upper and lower bounderies for optimization
constraintsLow = [eps 0 eps 0]; % [ N f sig rho]
constraintsUpp = [Inf 100 100 2*pi];

dt_STIM = process.stimulus.dt;
Fs_STIM = 1/dt_STIM;
time_STIM = linspace(0,  length(process.stimulus.Yuncorr)*dt_STIM,  size(process.stimulus.Yuncorr,1));
STA_LENGTH = 35;
process = process.CalcSTA( STA_LENGTH, 5);



%%
%         figure(1);
%         plot(time,targetfun(time));
%         figure(2);plot(time,modelfun(p(1),p(2),p(3),p(4)))
%         figure(1);plot(time,Y_est)
%         pause(1);
% %% Optimization 
% 
% % using the same optimization procedure as in section 2.1
% x0 = [1 1 1 1]; % start point
% time = process.STA.timeSTA;
% numcoeff = 4;
% N = length(time);
% dfe = N - numcoeff; % degrees of freedom - number of observations - num of coefficients
% STA = {process.STA.STA}';
% % modelfun = @(sig, f, mu, phi) exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ); % original Kernel model
% % costfun = @(Y) (@(x) sum( (Y' - modelfun(x(1),x(2),x(3),x(4))).^2 )); % SA with direct search  [sig f mu phi]
% modelfun = @(sig, f, mu, phi) normax( exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ) ); % with normax
% costfun = @(Y) (@(x) sum( ( Y' - modelfun(x(1),x(2),x(3),x(4))).^2 )); % SA with direct search  [sig f mu phi]
% 
% % [p_est, output, goodness, res, est] = deal(cell(process.numChannels,1));
% estSTA = cell(process.numChannels,1);
% %% Solving with a parfor Loop
% 
% [estSTA, res] = deal(zeros(process.numChannels, STA_LENGTH));
% % goodness = mat2cell(zeros(process.numChannels, 1));
% goodness = cellfun(@(x) struct('sse',[],'rsquare',0,'dfe',[],'adjrsquare',[],'rmse',[]), num2cell(1:process.numChannels));
% processArray = (1:process.numChannels); % the processes which needs to be estimated
% 
% % parameters for mulitple runs of the optimization
% AreAllOptim = true( process.numChannels, 1); % flag if the fitting has R2 > .85
% maxRun = 20; % maximum nubmbers of runs of the optimization
% numIterGlobal = 0; % index for iteration number of whole process
% 
% optimData = struct('LowConstr',LowConstr,'UppConstr',UppConstr,'numcoeff',...
%     numcoeff,'process',process,'STA_LENGTH',STA_LENGTH,'costfun',costfun,...
%     'modelfun',modelfun,'N',N,'dfe',dfe);
% [ estSTA, res, goodness, estPARAM ] = multiOptim_PARFOR( STA, optimData );
% 
% range = (1:process.numChannels);
% optimResults = struct('time',time,'data',STA,'estimationSTA',estSTA,...
%     'estimationPARAM',estPARAM,'Stat_goondness',goodness);               
% %% Parralel loop through all STAs
% % 2DO - each worker overrides 2 calculations (except of the last one)
% % generalize to the case where the local parrallel profile has less workers
% % than data
% 
% processArray = (1:process.numChannels); % the processes which needs to be estimated
% numIterD = distributed( processArray );
% % numIterD = distributed( (1:6) );
% 
% % parameters for mulitple runs of the optimization
% AreAllOptim = true( process.numChannels, 1); % flag if the fitting has R2 > .85
% maxRun = 3; % maximum nubmbers of runsof the optimization
% numIterGlobal = 0; % index for iteration number of whole process
% optimResults = cell(maxRun,1);
% 
% while any( AreAllOptim ) && numIterGlobal < maxRun % if R^2 > .85 and didn't exceed maxRun
%     spmd
%     %         [estSTA, res] = deal(codistributed.zeros(process.numChannels,STA_LENGTH));
%     %         goodness =  codistributed.cell(process.numChannels,1);
%             numIterLP = getLocalPart ( numIterD );
%             for n = numIterLP
%     %         for n = numIterLP(1):numIterLP(end)      
%                 fprintf('working on %d''th STA fitting ...\n', n);    
% 
%                 costIter = costfun(normax( STA{n} ));
%     %             [p_est,~,~,output] = annealData2Acorr(x0, LowConstr,UppConstr,costIter);          
%     %             [p_est,~,~,output] = gps_Data2STA(x0 , LowConstr,UppConstr, costIter);
%                 [p_est,~,~,output] = ga_Data2STA( numcoeff, LowConstr,UppConstr, costIter);
%     %             fprintf('current n is %d \n' ,n); 
%                 estSTA = modelfun( p_est(1), p_est(2), p_est(3), p_est(4)); % [sig f mu phi]   
%     %             plotFit( time, normax(STA{n}) , estSTA, 'STA 1 fit')
% 
%                 % Evaluating goodness of optimization and fit
%                 res = normax(STA{n}) - estSTA';  
%                 goodness = iGoodnessStructure(estSTA,[],res,dfe,N);
%                 
%                 if goodness.rsquare < .85
%                     AreAllOptim(n) = false;
%                 end
%     %             estSTA(n,:) = modelfun( p_est(1), p_est(2), p_est(3), p_est(4)); % [sig f mu phi]   
%     %             res(n,:) = normax(STA{n}) - estSTA';  
%     %             goodness(n) = iGoodnessStructure(estSTA,[],res,dfe,N);
%             end
%     end  
%     AreAllOptim = logical( prod(cell2mat(AreAllOptim(1:end)),2) ); % uniting the composite values to 
%     disp(['Are All Optimized? ',num2str( AreAllOptim' )])
%     numIterGlobal = numIterGlobal + 1;
%     numIterD = distributed(  processArray(AreAllOptim) );
%     optimResults{numIterGlobal} = struct('time',time,'data',{STA(range)},'estimation', {estSTA(range)},'Stat_goondness',{goodness(range)},'optimOut', {output(range)});               
% end
% range = (1:process.numChannels);
% optimResults = struct('time',time,'data',{STA(range)},'estimation', {estSTA(range)},'Stat_goondness',{goodness(range)},'optimOut', {output(range)});               
% 
% temp = goodness(range);cellfun(@(x)disp(x.rsquare),temp)
