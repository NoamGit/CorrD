function [ out ] = LoadDataIntoProcess( varargin )
%LoadDataIntoProcess() this function gets data files, and its type and build a data
%structure dedicated for Neural singnal Analysis
%   Detailed explanation goes here

% handle different kind of inputs
switch nargin
    case 0 
        [ filename, pathname ] = uigetfile('*.mat', 'Select m-fil data');
        currentFolder = pwd;
        cd(pathname);
        varargin = importdata(filename);
        cd(currentFolder); 
       
        % 2DO : here should come a header read but I only have Ronens data for
        % know
        signalArray = varargin;
        type = 'fileload';
        
    otherwise % default data set for Segev's data
        % recall that in practice we are running 15 trails of 60 sec
        % experiment on all 23 cells. This can yield way more data. Further
        % D.W(:,1) probabily indicates on times that the stimulus is "refreshed"
        % so it might be a an empty gap 
        
%         D = load('C:\Users\Noam\OneDrive - Technion\Project - CorrD\Code and toolboxes\Data\SpikeTimeGauss1B.mat');
        D = load('C:\Users\noambox\Documents\CorrD\SourceData\Ronen\SpikeTimeGauss1B.mat'); % change path
        stimRate = 30; % 30 [hz]        
        sr = 10e3; % sampling rate [hz]
        sr_new = 450; % resampling parameter [hz] rem( 480 , 30) = 0 ;
        modelfunc = @(mu, sig) exp( mu + sig * t);
        exp_model = struct('type', 'exp', 'model', modelfunc,'estimation',[] ); % default kernel is exp
        spiketimes_all = cellfun(@(x) x .* 1/sr, {D.TT.sp},'UniformOutput',false);
        spiketimes = cellfun(@(x) x( x <= 900 ), spiketimes_all,'UniformOutput',false); % only first 900 sec are used for this experiment
        stimulusRAW = D.W(:,3); % W(:,3) matrix
        
        % take 15 trails of a 60 sec experiment (900 sec for each cell)
        sTHRESCORR = 15 * 60 * stimRate; % samples [ numTrails * lengthOfExperiment * stimulusRate ]
        stimulusUC = stimulusRAW(1:sTHRESCORR);
        stimulusUC = stimulusUC - mean(stimulusUC); % for uncorrelated WGN stimulus
        
        % build object
        signalArray = PointProcessData( sr, spiketimes,stimRate, stimulusRAW, sr_new, exp_model, stimulusUC);
        type = 'Segev';
end
        out = struct('type', type, 'Data',signalArray);
end
