classdef PointProcess
%%======================================= Linear kernel model estimation===========================================
% In this class we simulate the retrival of a Linear Kernel.
% We go back and fort with the simulation and see if the estimation yields good results
% =====================================================================================================
       
    properties(GetAccess = 'public', SetAccess = 'public')
        % Raw properties
        Fs  % sample rate [Hz]
        dt  % time resolution [s]
        N   % num of sampels
        T   % signal length [s]
        t   % timeline 
        
        % Statistical properties
        acorr % autocorrelation functions
        
        % model properties
        CGP % correlated Gaussian Process
        lambda % Rate Process
        spiketimes % spiketimes of the point process
        linKernel % linear Kernel model properties
        NlinKernel % Non Linearity
        stimulus % stimulus
    end
    
    methods
        % constructor
        function obj = PointProcess(Fs, N)
            if nargin > 0
                obj.Fs = Fs;   
                obj.dt = 1/Fs;  
                obj.N = N;  
                obj.T = N/Fs;   
                obj.t = linspace(0,obj.T,N);
            end
        end
        
        function obj = set.linKernel( obj, model )
            obj.linKernel = model;
        end 
        
        function obj = set.NlinKernel( obj, model )
            obj.NlinKernel = model;
        end
        
        [mfr] = ppmean(obj)
    end
    
end

