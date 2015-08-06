classdef PointProcessSynt < PointProcess
%%=======================================  PointProcessSynt ===========================================
% In this class we simulate a point process in a Forward scheme. This uses
% the LNP model, and we can estimate the STA acorr and other statistical
% properties of the signal
% We go back and fort with the simulation and see if the estimation yields good results
% =====================================================================================================
        
    properties
        % NLP model properties
        noiseMean % parameters of the noise
        noiseStd
        
        % Parameter estimation
        PEresults % results of parameter estimation
    end
    
    methods
        function obj = PointProcessSynt(Fs, N, noiseMean, noiseSTD)
            if nargin > 0
                obj.Fs = Fs;   
                obj.dt = 1/Fs;  
                obj.N = N;  
                obj.T = N/Fs;   
                obj.t = linspace(0,obj.T,N);
                obj.noiseMean = noiseMean;      
                obj.noiseStd = noiseSTD;
                obj.stimulus = noiseMean + noiseSTD * randn(N, 1);
            end
        end
                               
        function obj = generateProcess( obj, amplitude )
            % generateProcess( obj, amplitude ) generates process according to NLP
            timeSupport = obj.t(obj.linKernel.support);
            obj.CGP = conv(obj.stimulus , amplitude * obj.linKernel.model( timeSupport ),'same');
            obj.lambda = obj.NlinKernel.model( obj.CGP ); 
            obj.spiketimes = simpp( obj.lambda , obj.dt);  
            
            % generate statistical properties
            acorrCGP = xcorr( obj.CGP, obj.linKernel.maxlags,'unbiased' );
            acorrLambda = xcorr( obj.lambda, obj.linKernel.maxlags,'unbiased');
            acorrPP =  ppcorr(obj.spiketimes, obj.dt, obj.linKernel.maxlags, 'unbiased', 0);
            obj.acorr = struct('CGP_Corr', acorrCGP,'lambdaCorr', acorrLambda,'PointP_Corr', acorrPP);
        end
        
        function obj = estimatesStatsProperties( obj, maxlags, varargin )
            % given     estimates the acorr of the estimated CGP
            % if value in centerLag is given we extract it from varargin -
            % the acorr for lambda
            if isempty(obj.spiketimes)
                error('No spiketimes where loaded');
            else
                obj.acorr.lambdaCorr =  ppcorr(obj.spiketimes, obj.dt, maxlags, 'unbiased', 1); % doesn't take value in lag 0
                if nargin > 2 % execute if zeroLag is given 
                    zeroLag = ceil( length(varargin{1})/2 );
                    obj.acorr.lambdaCorr(zeroLag) = varargin{1}(zeroLag); 
                else
                    zeroLag = ceil( length(obj.acorr.lambdaCorr)/2 );
                end
                E_lambda = ppmean( obj.spiketimes );
                % See Krumin's cacls in p.20
                sigma = sqrt(log(obj.acorr.lambdaCorr(zeroLag)/(E_lambda.^2)));
                obj.acorr.CGP_Corr = (1/sigma^2)*log(obj.acorr.lambdaCorr/(E_lambda.^2));    
            end
        end
        
        [] = plotEstimation( obj, data )
        
        obj = paramEstim( obj, Data, time, Y_hat, varargin )
         
    end
    
end

