classdef ComplexCellSimulation < PointProcess
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % stimulus parameters
        noiseMean % parameters of the noise
        noiseStd
        
        X
        nl
        lambd_Uni
        CP 
        
        filterMatrix
        numOfFilters
        filterRelationfunc
    end
    
    methods
        function obj = ComplexCellSimulation( varargin )
             % initializes stimulus
        if any(isstruct( varargin{:} ))
            if( any(isfield(varargin{1},'Fs')) );obj.Fs = varargin{1}.Fs; end;
            if( any(isfield(varargin{1},'size')) );obj.N = varargin{1}.size; end;
            if( any(isfield(varargin{1},'noiseMean')) );obj.noiseMean = varargin{1}.noiseMean; end;
            if( any(isfield(varargin{1},'noiseSTD')) );obj.noiseStd = varargin{1}.noiseSTD; end;
            if( any(isfield(varargin{1},'nl')) );obj.nl = varargin{1}.nl; end;
            if( any(isfield(varargin{1},'filterMatrix')) );obj.filterMatrix = varargin{1}.filterMatrix; end;
            if( any(isfield(varargin{1},'filterRelationfunc')) ); obj.filterRelationfunc = varargin{1}.filterRelationfunc; end;
        end 
            obj.dt = 1/obj.Fs;  
            obj.numOfFilters = size(obj.filterMatrix,2);
            obj.T = obj.N/obj.Fs;   
            obj.t = linspace(0,obj.T,obj.N);
            obj.stimulus = obj.noiseMean + obj.noiseStd * randn(obj.N, 1);
         end
                               
        function obj = generateProcess( obj, amplitude )
            % generateProcess( obj, amplitude ) generates process
            [ obj.X, obj.lambda ] = deal( zeros( length( obj.stimulus ), obj.numOfFilters) );
            for k = 1:obj.numOfFilters
                cgp = conv( obj.stimulus, amplitude * flipud( obj.filterMatrix(:,k) ), 'full');
                obj.X(:,k) = cgp(1:length( obj.stimulus ));
                obj.lambda(:,k) = obj.nl(obj.X(:,k));
            end
            
            obj.lambd_Uni = obj.filterRelationfunc(obj.lambda); 
            obj.spiketimes = simpp( obj.lambd_Uni , obj.dt);
            obj.CP = histc(obj.spiketimes, obj.t);
            disp(['simulated mean spike rate - ', num2str(numel(obj.spiketimes)/( length( obj.stimulus ) * obj.dt ) ), ' [Hz]'])
        end
        
        
    end
    
end

