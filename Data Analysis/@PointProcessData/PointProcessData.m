classdef PointProcessData < PointProcess
    %PointProcessData a point process class that we load data into (not simulate)
    %   Detailed explanation goes here
    
    properties
        numChannels
        maxlags
        
        % features
        CP % Couting process from spiketimes
        STA % real and estimated STA
    end
    
    methods
        %% default constructor
        function obj = PointProcessData(Fs, spiketimes, stimRate ,stimulus, varargin)
            % in varargin you can:          change sampling rate
            %                               change nonlinear filter
            if nargin > 0
                
                % Preprocessing (1) Resampling - computatinally and biologically working with oversampled
                % data is difficult and inefficient. 
                % TODO: stimulus and the respone to have the same rate.
                if nargin > 4
                    newrate = varargin{1};
                    obj.NlinKernel = varargin{2};
                    stimulusUC = varargin{3};
                else 
%                     newrate = Fs;
%                     exp_model = @(mu, sig) exp( mu + sig * t);
%                     obj.NlinKernel = struct('type', 'exp', 'model', value2 ); % default kernel is exp
                end
                [ Fs_out,spiketimes_out ] = resampleSpikeTimes( Fs, spiketimes, newrate );
                
                obj.Fs = Fs_out; 
                obj.dt = 1/Fs_out;  
                obj.spiketimes = spiketimes_out; % is a cell array
                obj.stimulus = struct('Yraw', stimulus,'Yuncorr',stimulusUC,'dt', 1/stimRate);
                obj.T = length(stimulusUC)*1/stimRate; % end of stimulus
                obj.N = ceil( obj.T * obj.Fs );  
                obj.t = linspace(0,obj.T,obj.N);
                obj.CP = cellfun(@(x) histc(x, obj.t), obj.spiketimes,'UniformOutput',false);
                obj.numChannels = numel(spiketimes);
                obj.acorr = struct('lags',[],'CGP_Corr_RAW', [],'lambda_Corr_RAW', [],'CGP_Corr_FILT', [],'lambda_Corr_FILT', []);
            end
        end
        
        %% Calculate STA of process
        function obj = CalcSTA( obj, varargin )
            % uses compute_sta to find the Spike triggered averge of the PP
            % out of the stimulus and the CP. STA duration should be
            % provided. We also upsample the stimulus in case it has
            % different SR from the CP.
            %   varg{1} STA_NUMSAMPLES - STA duration
            %   varg{2} STA_LWE - STA left window edge
            %   varg{3} resampleFlag - flag for resampling the stimulus to the sp
            %                   resolution
            % ** note: run first PreProcessCP to obtain the CP
            
            % check exsistance of data 
            assert(~isempty(obj.stimulus), '     no CP or stimulus loaded for the process')
            
            resampleFlag = any(strcmp(varargin,'resample stimulus'));
            numArg = nargin-1;
            
            if numArg < 2
                STA_NUMSAMPLES = 35 ; % [ms] default STA of 35 samples 
                STA_LWE = 5;
            else
                STA_NUMSAMPLES = varargin{1};
                STA_LWE = varargin{2};
            end
            
            % resample if flagged
            if resampleFlag
                upsamplingFactor = obj.Fs/(1/obj.stimulus.dt); % must be  
                assert( ~rem(obj.Fs,(1/obj.stimulus.dt)) ,...
                    '  resampling failed. The F of stim and process have non zero remainder');
                
                obj.stimulus.Yuncorr = obj.stimulus.Yuncorr(1:ceil(obj.t(end)/obj.stimulus.dt)); % prunning only to stimulus in the time axis
                upsampledStim = repmat(obj.stimulus.Yuncorr',[upsamplingFactor 1]);
                obj.stimulus.Yuncorr = upsampledStim(:);
                obj.stimulus.dt = obj.dt;
                STA_NUMSAMPLES = STA_NUMSAMPLES * upsamplingFactor;
                STA_LWE = STA_LWE * upsamplingFactor;
            end
            
            % use compute sta func to evalate STA
            if iscell(obj.spiketimes) % multi process
                [ sta ] = cellfun(@(ST) compute_sta( obj.stimulus.Yuncorr, ST, STA_NUMSAMPLES, STA_LWE, obj.stimulus.dt)...
                    ,obj.spiketimes,'UniformOutput',false);
            else % single process
                [ sta ] = compute_sta( obj.stimulus.Yuncorr, obj.spiketimes, STA_NUMSAMPLES, STA_LWE, obj.stimulus.dt );
            end
            STA_TIME = (-(STA_NUMSAMPLES-STA_LWE-1):STA_LWE) * obj.stimulus.dt;
            obj.STA = struct('timeSTA',STA_TIME,'realSTA',sta,'estSTA',[]);
        end
        
        %% proprocess data according to LNP model
        obj = PreProcessCP( obj, varargin )
        
        %% find PP correlation
        [lambda_corr] = ppcorr(obj, bias, correctR0, method)
        
        %% estimate non linearity
        obj = nlestimation( obj, varargin );
        
        %% show spike time interval histogram and exponential fit
        showISI( obj, varargin ); 
        
    end
    
end

