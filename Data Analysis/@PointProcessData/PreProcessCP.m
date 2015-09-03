function obj = PreProcessCP( obj,varargin ) 
%PreProcessData() here we restore from spiketimes the following features:
% the CGP, CGP acorrRaw and CP acorrSmooth
%       inputs:     flagCompare = varargin{1};
%                   cutoff freq = varargin{2}; 

% [ CGP_Corr_RAW, CGP_Corr_FILT, lambda_Corr_FILT, lambda_Corr_RAW ] = deal(cell(obj.numChannels,1));

% define the functions for moving to CGPs
flagCompare = 0;
if nargin > 1 %define the non linearity type
    NLtype = varargin{1};
    flagDoublet = varargin{2};
    flagCompare = varargin{3};
else % default settings
    NLtype = 'exp';
end

% ** optional - Doublet cancelation for Poissonian correction
if flagDoublet
    obj = DoubletCancel(obj);
end

% mean firing rate calculation using Krumin's function 
% (We can change to calculate all commands only on one side of the symmetrical AC) 
E_lambda = ppmean(obj.spiketimes); 
lambda_Corr_RAW = ppcorr( obj, 'unbiased', 1, 'auto' ); % moment generating function for estimating value in tau=0 
EStack = cell2mat(E_lambda);
corrStack_L = cell2mat(lambda_Corr_RAW') ;

% calculate CGP for modeldriven non linearity
if( any(strcmp(NLtype,{'abs','sqrt','exp'})) )
    CGPestimation_ModelDriven( obj, NLtype )
end

% Assign data to Data structure
obj.acorr.lags = (-obj.maxlags:obj.maxlags) .* obj.dt;
obj.acorr.lambda_Corr_RAW = mat2cell( corrStack_L, length(corrStack_L) , ones(obj.numChannels,1)' );

% compare visually to extract good data
if flagCompare
    lambdaCorr_forw = cellfun(@(x1,x2) xcorr( polyval(x1,x2), obj.maxlags,'unbiased') ,...
                        {obj.NlinKernel.estimation(:).polyfit}', obj.CGP(:),...
                        'UniformOutput',false);                    
    lambdaCorr_back = obj.acorr.lambda_Corr_RAW;   
    plotCompare( lambdaCorr_back, lambdaCorr_forw, obj.acorr.lags, 6, 'cell', 1, (1:24));
end

end


