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
       
        % FIX : here should come a header read but I only have Ronens data for
        % know
        signalArray = varargin;
        type = 'fileload';
        
    otherwise % default data set for Segev's data
        D = load('C:\Users\Noam\OneDrive - Technion\Project - CorrD\Code and toolboxes\Data\SpikeTimeGauss1B.mat');
%         D = load('D:\# Projects (Noam)\# CorrD\Code and toolboxes\Data\SpikeTimeGauss1B.mat');
        stimRate = 30; % 30 [hz]        
        sr = 10e3; % sampling rate [hz]
        sr_new = 500; % resampling parameter [hz]
        modelfunc = @(mu, sig) exp( mu + sig * t);
        exp_model = struct('type', 'exp', 'model', modelfunc ); % default kernel is exp
        spiketimes = cellfun(@(x) x .* 1/sr, {D.TT.sp},'UniformOutput',false);
        stimulusRAW = D.W(:,3); % W(:,3) matrix
        
        % uncorrelating Segev's threshold
        sTHRESCORR = 5.346e4; % stimulus threshold correlation
        stimulusUC = stimulusRAW(1:sTHRESCORR);
        stimulusUC = stimulusUC - mean(stimulusUC); % uncorrelated WGN stimulus

        signalArray = PointProcessData( sr, spiketimes,stimRate, stimulusRAW, sr_new, exp_model, stimulusUC);
        type = 'Segev';
end
        out = struct('type', type, 'Data',signalArray);
end

