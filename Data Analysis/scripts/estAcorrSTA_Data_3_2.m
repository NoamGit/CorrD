% estAcorrSTA_Data_3_2 

%% load data module
% comparing parametric model to fit data in forward and backward
% procedure. The comparison is done in the nl and the CGP plane

% Loads data (also calculates CP)
    nlType = 'exp';
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
    flagDoublet = 0;
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
    showFlag = 0;
    process = process.nlestimation( poly_order, showFlag, amp, nlType,'original stimulus' ); 

% ** optional - see function description - result : exlude cells 9 
    compareInNlPlane( process );

% ** optional - see func description
% result - good statistical prop in cell 1 2 4 5 7 8 9 14 15 20 21 22  
    showISI( process );

% ** optional - 
    burstDetection( process );
    
%% STC analysis module


[ sta, w_sta ] = compute_white_sta(process.stimulus, process.CP, 3/(1/Fs) );

[e_vec, est_sta, var_prec] = compute_stc( cpxproc.stimulus, cpxproc.CP, 3/(1/Fs), sta, 'suppress sta');
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



