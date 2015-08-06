function obj = paramEstim( obj, Data, t, varargin )
        % model based estimation of parameters from the CGP acorr
        % function. 
        % inputs:       Y_hat - hypotized model
        %               varargin{1} - lower optimization boundaries [ N f sig ]
        %               varargin{2} - upper optimization boundaries [ N f sig ]
        %               varargin{3} - FIR filter's coefficients for noise removal before estimation
                        
    if isempty(obj.acorr.CGP_Corr)
            error('No data for estimation');
    else
        
        % filters data if coefficients are provided
        if nargin > 3
            Y = filtfilt(varargin{3}, 1, obj.acorr.CGP_Corr); %  filter data with provided FIR filter
        else 
            Y = obj.acorr.CGP_Corr; % unfilterd Data 
        end
        
        % set model equation  
        fun = @(N,f,s) N * (exp((-t.^2)/(2*s^2)).*cos((2*pi*f).*t))-exp((-t.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0);
        costfun = @(p) sum( (Y' - fun(p(1),p(2),p(3))).^2 ); % SA with direct search
        if nargin > 2 
            LowConstr = varargin{1}; % low optimization boundaries [ N f sig ]
            UppConstr = varargin{2}; % upp optimization boundaries
        else
            LowConstr = [-Inf 0 0]; % low optimization boundaries
            UppConstr = [Inf 40 1e3]; % upp optimization boundaries
        end
        
        % uses matlab's Simmulated annealing with pattern search
        % otimization
        paramData = Data.linKernel.parameters( [1,4] );
        x0 = [1 1 1]; % start point
        fprintf('starting optimization with SA ...\n');
        [p,~,~,output] = annealData2Acorr(x0, LowConstr,UppConstr,costfun);
        Y_est = fun( p(1), p(2), p(3));
        numcoeff = 3;
        N = length(t);
        dfe = N - numcoeff; % degrees of freedom - number of observations - num of coefficients
        res = Y' - Y_est;
        goodness = iGoodnessStructure(Y_est,[],res,dfe,N);
        optimResults = struct('estimation', Y_est,'optimOut', output);
        
        if strcmp( Data.linKernel.type , 'gaussSine' )
            paramEst = abs([p(3),p(2)]); % [sig f]
            paramData = Data.linKernel.parameters( [1,4] );
            param = (paramData - paramEst)./mean([paramData;paramEst]);
            est_RMS = sqrt( (param * param' )/numel(param) ); % RMS of the estimation
        else
            est_RMS = [];
        end
        
        obj.PEresults = struct('OptimResults',optimResults,'GoodnessStatistics',goodness,'paramData',paramData,...
            'paramEst',paramEst,'estimationRMS', est_RMS);
    end
end

