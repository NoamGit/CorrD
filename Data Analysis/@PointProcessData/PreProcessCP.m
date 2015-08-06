function obj = PreProcessCP( obj,varargin ) 
%PreProcessData() Here we take the spiketimes and the general data structure and compute the following features:
% the CP, CP acorrRaw and CP acorrSmooth
%       inputs:     NLtype = varargin{1};
%                   cutoff freq = varargin{2}; 

[ CGP_Corr_RAW, CGP_Corr_FILT, lambda_Corr_FILT, lambda_Corr_RAW ] = deal(cell(obj.numChannels,1));

% define the functions for moving to CGPs
if nargin > 1 %define the non linearity type
    NLtype = varargin{1};
    cutoff = varargin{2};
else % default settings
    NLtype = obj.NlinKernel.type;
    cutoff = 8; % cutoff freq hz
end

switch NLtype
        case 'exp' % (Exponential Non Linearity)
            % moving to CGP data 
            findSigma = @(X,Y) sqrt( log(X(obj.maxlags+1,:) ./ Y.^2) ) ; % X is the lambda_AC and Y the mean FR (p.26)
            findMu = @(X,Y) log( Y.^2./ X(obj.maxlags+1) ) ;
            findCGP = @(X,Y) bsxfun(@times ,log( bsxfun(@rdivide , X, Y.^2) ),( 1./findSigma(X,Y).^2 ) ); % p.26
%               Int_AC{n}(abs(Int_AC{n}) < 0.001) = min(Int_AC{n}(Int_AC{n} > 1e-5)) * 0.01; % By restraining Int_AC to a number we are avoiding the log transform to a complex value
%               sigma = (log(Int_AC{n}(idx_tau_0)/E_lambda{n}.^2)); % estimated sigma 
%               CGP{n} = (1/sigma.^2)*log(Int_AC{n}/E_lambda{n}.^2); % Correlation pre Distortionaccording to michael cacls in p.20  

        case 'sqr' % (Square value Non Linearity) FILL ACCORDING TO p.26
    
%              [mu sigma] = msfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}),'square'); % michaels function for abs NL (inspired by p. 39-40)
%              CGP{n} = grfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}), mu, sigma,'square'); 

        case 'abs' % (Absolute Non Linearity)
%               Int_AC{n}(abs(Int_AC{n}) < 0.001) = min(Int_AC{n}(Int_AC{n} > 1e-5)) * 0.01; % By restraining Int_AC to a number we are avoiding the log transform to a complex value
%               [mu sigma] = msfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}),'abs'); % michaels function for abs NL (inspired by p. 39-40)
%               CGP{n} = grfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}), mu, sigma,'abs');
end

% mean firing rate calculation using Krumin's function (We can change to calculate all commands only on one side of the symmetrical AC) 
E_lambda = ppmean(obj.spiketimes); 
lambda_Corr_RAW = ppcorr( obj, 'unbiased', 0, 'auto'); % here we use the moment generating function for estimating value in tau=0 
EStack = cell2mat(E_lambda);
corrStack_L = cell2mat(lambda_Corr_RAW');
faultyVal = corrStack_L <=0;
corrStack_L( faultyVal ) = min(corrStack_L(corrStack_L>0));
corrStack_CGP = findCGP( corrStack_L, EStack );
sigSize = obj.maxlags*2+1;
fv_linInd = find(faultyVal);
fv_edge = mod( fv_linInd, sigSize ) == 1 | mod( fv_linInd, sigSize ) == 0;
[t1] = interpolateFV(fv_linInd(~fv_edge)', corrStack_CGP);
% corrStack_CGP(fv_linInd(~fv_edge)) = 


for ii = 1:obj.numChannels
    
    % assure definite positive signal
    
    
    % compute the CGP Corr
%     CGP_Corr_RAW{ii} = arrayfun findCGP( lambda_Corr_RAW{ii}, E_lambda{ii} );
    
    % filter high frequencies - FIX
    lambda_Corr_FILT{ii} = LowPassAcorr( obj, lambda_Corr_RAW{ii}, cutoff );
    CGP_Corr_FILT{ii} = LowPassAcorr( obj, CGP_Corr_RAW{ii}, cutoff );
end

    % assign
%     obj.acorr.
end

    function [ out ] =interpolateFV( ind, values )
    % function which interpolates faulty alab values from AC calculation by
    % takung the mean of their closure
    
        q = ind;
        
        % find consequtive indexs in array
        q_tag = [0 q inf];
        a=diff(q_tag);
        b=find( [a inf]>1 );
        c=diff([0 b]); % length of the sequences
        d=cumsum(c); % endpoints of the sequences
        
        d_p = [1 d(1:end-1)+1];
        IND = d-d_p ~= 0;
        fvClosureMeanV = arrayfun( @(x,y) 0.5 * (values(x-1)+values(y+1)), q_tag(d_p(IND)), q_tag(d(IND)) );
        out = [ fvClosureMeanV' q_tag(d_p(IND))' q_tag(d(IND))' ];
    end
