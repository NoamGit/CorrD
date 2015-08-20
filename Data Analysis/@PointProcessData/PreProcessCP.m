function obj = PreProcessCP( obj,varargin ) 
%PreProcessData() here we restore from spiketimes the following features:
% the CGP, CGP acorrRaw and CP acorrSmooth
%       inputs:     flagCompare = varargin{1};
%                   cutoff freq = varargin{2}; 

temp = 0;
[ CGP_Corr_RAW, CGP_Corr_FILT, lambda_Corr_FILT, lambda_Corr_RAW ] = deal(cell(obj.numChannels,1));

% define the functions for moving to CGPs
flagCompare = 0;
if nargin > 1 %define the non linearity type
    NLtype = varargin{1};
    flagCompare = varargin{2};
else % default settings
    NLtype = exp;
end

switch NLtype
        case 'exp' % (Exponential Non Linearity)
            % moving to CGP data 
            findSigma = @(X,Y) sqrt( log(X(obj.maxlags+1,:) ./ Y.^2) ) ; % X is the lambda_AC and Y the mean FR (p.26)
            findMu = @(X,Y) log( Y.^2./ X(obj.maxlags+1) ) ;
            findAcorrCGP = @(X,Y) bsxfun(@times ,log( bsxfun(@rdivide , X, Y.^2) ),( 1./findSigma(X,Y).^2 ) ); % p.26
            typeFlag = 1;

        case 'sqr' % (Square value Non Linearity) FILL ACCORDING TO p.26
    
%              [mu sigma] = msfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}),'square'); % michaels function for abs NL (inspired by p. 39-40)
%              CGP{n} = grfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}), mu, sigma,'square'); 
%              typeFlag = 1;

        case 'abs' % (Absolute Non Linearity)
%               Int_AC{n}(abs(Int_AC{n}) < 0.001) = min(Int_AC{n}(Int_AC{n} > 1e-5)) * 0.01; % By restraining Int_AC to a number we are avoiding the log transform to a complex value
%               [mu sigma] = msfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}),'abs'); % michaels function for abs NL (inspired by p. 39-40)
%               CGP{n} = grfind(num2cell(E_lambda{n}), num2cell(Int_AC{n}), mu, sigma,'abs');
%               typeFlag = 1;    

        otherwise % estimate non linearity from data
            typeFlag =0;
%           findSigma = @(X,Y) sqrt( log(X(obj.maxlags+1,:) ./ Y.^2) ) ; % X is the lambda_AC and Y the mean FR (p.26)
%           findMu = @(X,Y) log( Y.^2./ X(obj.maxlags+1) ) ;
%           findCGP = @(X,Y) bsxfun(@times ,log( bsxfun(@rdivide , bsxfun(@minus, X ,-1*ones(1,23)), Y.^2) ),( 1./findSigma(X,Y).^2 ) ); % p.26
end

% mean firing rate calculation using Krumin's function 
% (We can change to calculate all commands only on one side of the symmetrical AC) 
E_lambda = ppmean(obj.spiketimes); 
lambda_Corr_RAW = ppcorr( obj, 'unbiased', 0, 'auto'); % moment generating function for estimating value in tau=0 
EStack = cell2mat(E_lambda);
corrStack_L = cell2mat(lambda_Corr_RAW') ;

% fix and interpolate numerical faulty values of the AC calculation
if(typeFlag)
    faultyVal = corrStack_L <= 1e-5; % faulty values are zero approaching
    corrStack_L( faultyVal ) = min(corrStack_L(corrStack_L>0));
    corrStack_CGP = findAcorrCGP( corrStack_L, EStack );
    sigSize = obj.maxlags*2+1;
    fv_linInd = find(faultyVal); % faulty vals linear index
    fv_edge = mod( fv_linInd, sigSize ) == 1 | mod( fv_linInd, sigSize ) == 0;
    fv_start  = fv_linInd( mod( fv_linInd, sigSize ) == 1 );
    fv_end  = fv_linInd( mod( fv_linInd, sigSize ) == 0 );
    corrStack_CGP( fv_start ) = 0; % fix all edge fv to 0
    corrStack_CGP( fv_end ) = 0; % fix all edge fv to 0
    interpolateFV(fv_linInd(~fv_edge)); % interpolant all non edges faulty values
    
    % Apply smoothing & assign value
    obj.acorr.CGP_Corr_RAW = mat2cell( corrStack_CGP, length(corrStack_CGP) , ones(obj.numChannels,1)' );
    obj.acorr.CGP_Corr_FILT = cellfun( @(X) smooth(X, 'sgolay',4), CGP_Corr_RAW, 'UniformOutput' ,false);
end
obj.acorr.lags = (-obj.maxlags:obj.maxlags) .* obj.dt;
obj.acorr.lambda_Corr_RAW = mat2cell( corrStack_L, length(corrStack_L) , ones(obj.numChannels,1)' );

% compare visually to extract good data
if flagCompare
    lambdaCorr_forw = cellfun(@(x1,x2) xcorr( polyval(x1,x2), obj.maxlags,'unbiased') ,...
                        {obj.NlinKernel.estimation(:).polyfit}', obj.CGP(:),...
                        'UniformOutput',false);
                    
    lambdaCorr_back = obj.acorr.lambda_Corr_RAW;   
    % What if the correlated noise effect was only convolving with the
    % g = Yuncorr_correlation
%     g = xcorr(obj.stimulus.Yuncorr,obj.maxlags,'unbiased');
%     G = xcorr(g,'unbiased');
%     lambdaCorr_back = cellfun(@(x) conv( x, G,'same'),obj.acorr.lambda_Corr_RAW,'UniformOutput',false);
    plotCompare( lambdaCorr_back, lambdaCorr_forw, obj.acorr.lags, 6, 'cell', 1, (1:24));
end

    function [ ] = interpolateFV( ind )
    % function which interpolates faulty acorr values from AC calculation by
    % taking the mean of their neighb
        numIdx = length(ind);
        
        % extract lenght of conseq series
        ind_pad = [0 ;ind; inf];
        a=diff(ind_pad);
        b = find( [a ;inf]>1 );
        b_tag = [a; inf] > 1;
        b_tag = b_tag & circshift(b_tag, -1); % all non conseq fv indexs
        c = diff([0 ;b]); % all non 1 are sequence
        c_idx = c ~= 1; % location with the length of the sequence  
        d=cumsum(c)-1; % endpoints of the sequences corresponding to c
        
        % interpolate single fv indecators and fv series
        corrStack_CGP(ind(b_tag(1:numIdx))) = mean([corrStack_CGP(ind(b_tag(1:numIdx))-1);...
            corrStack_CGP(ind(b_tag(1:numIdx))+1)]); % means non conseq fv indexs
        for  k = find(c_idx)' % fv conseq series
            corrStack_CGP(ind(d(k))+1-c(k):ind(d(k))) = mean(corrStack_CGP([ind(d(k))-c(k), ind(d(k))+1]));
        end
    end 
end


