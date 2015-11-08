function [ out ] = LoadDataIntoProcess( varargin )
%LoadDataIntoProcess() this function gets data files, and its type and build a data
%structure dedicated for Neural singnal Analysis
%   Detailed explanation goes here

% handle different kind of inputs
switch varargin{1}
    case '' 
        [ filename, pathname ] = uigetfile('*.mat', 'Select m-fil data');
        currentFolder = pwd;
        cd(pathname);
        varargin = importdata(filename);
        cd(currentFolder); 
       
        % 2DO : here should come a header read but I only have Ronens data for
        % know
        signalArray = varargin;
        type = 'fileload'; 
        out = struct('type', type, 'Data',signalArray);
        
    case 'synt'
        % load synthesized data
        [ sr, stimRate ] = deal(1/varargin{2});      
        modelfunc = varargin{3}(1).nonlinKernel.model;
        nl_DATA = {varargin{3}.nonlinKernel};
        nl_model = struct('type', '2 degree poly',...
            'model', modelfunc,'estimation',[],'original', nl_DATA); % default kernel is exp
        spiketimes = {varargin{3}.spiketimes};
        stimulus = varargin{4};
        
        % build object
        out = PointProcessData( sr, spiketimes,stimRate, [], [], nl_model, stimulus);
        
%     case 'synt'
%         % self generated data
%         dt_stim = varargin{2}; 
%         samplingRate = 1/dt_stim; % [hz]
%         modelfunc = cellfun(@(x) x.model,{varargin{3}.nonlinKernel},'UniformOutput',false);
%         model_nl = struct('type', 'estimated', 'model', modelfunc,'estimation',[] ); % default kernel is exp
%         spiketimes = {varargin{3}.spiketimes}; % only first 900 sec are used for this experiment
%         stimulusUC = varargin{4}; % W(:,3) matrix
%         
%         % build object
%         signalArray = PointProcessData( samplingRate, spiketimes,samplingRate, [], [], model_nl, stimulusUC);
%         out = signalArray;
        
    otherwise % default data set for Segev's data
        % recall that in practice we are running 15 trails of 60 sec
        % experiment on all 23 cells. This can yield way more data. Further
        % D.W(:,1) probabily indicates on times that the stimulus is "refreshed"
        % so it might be a an empty gap 
        
        try
            D = load('C:\Users\noambox\Documents\CorrD\SourceData\Ronen\SpikeTimeGauss1B.mat'); % change path  
        catch exception
            try
                 D = load('C:\Users\Noam\Documents\GitHub\CorrD\SourceData\Ronen\SpikeTimeGauss1B.mat'); % change path
            catch exception
                 disp(exception.message);
                display('file was not found. Please load manually...');
                [fileN, pathN ] = uigetfile('*.mat','Select Data file');    % C:\Users\Noam\Documents\GitHub\...SpikeTimeGauss1B.mat
                D = load([ pathN, fileN ]); % change path
            end
        end        
        stimRate = 30; % 30 [hz]        
        sr = 10e3; % sampling rate [hz]
        sr_new = 450; % resampling parameter [hz] rem( 480 , 30) = 0 ;
        
        if nargin > 1
            switch varargin{2}
                case 'exp'
                    modelfunc = @(theta) theta(2) * exp( theta(1) .* t );
                    nl_model = struct('type', 'exp', 'model', modelfunc,'estimation',[] ); % default kernel is exp
                case 'poly'
                    modelfunc = @(theta) polyval(fliplr(theta), t);
                    nl_model = struct('type', 'poly', 'model', modelfunc,'estimation',[] ); % default kernel is exp
            end
        end
        
        spiketimes_all = cellfun(@(x) x .* 1/sr, {D.TT.sp},'UniformOutput',false);
        spiketimes = cellfun(@(x) x( x <= 900 ), spiketimes_all,'UniformOutput',false); % only first 900 sec are used for this experiment
        stimulusRAW = D.W(:,3); % W(:,3) matrix
        
        % take 15 trails of a 60 sec experiment (900 sec for each cell)
        sTHRESCORR = 15 * 60 * stimRate; % samples [ numTrails * lengthOfExperiment * stimulusRate ]
        stimulusUC = stimulusRAW(1:sTHRESCORR);
        stimulusUC = stimulusUC - mean(stimulusUC); % for uncorrelated WGN stimulus
        
        % build object
        signalArray = PointProcessData( sr, spiketimes,stimRate, stimulusRAW, sr_new, nl_model, stimulusUC);
        type = 'Segev';
        out = struct('type', type, 'Data',signalArray);
end
end

