function [ est ] = estSTA_Sim_2_2
% [est, proc] = estSTA_Sim_2_2 This "function" simulates multiple values for the linear Kernel and checks the
% goodness of fit to a predefined model 

% define parameter space to draw from
space_sig = linspace(5, 10 , 200);       space_mu = linspace(0, 10, 200); 
space_phi = linspace(0, pi, 200);         space_f = linspace(-0.5, 0.5 ,200); % [rad] & [hz] range abs[0.005, 2*Fs]

% define upper and lower bounderies for optimization
LowConstr = [eps space_f(1) space_mu(1) space_phi(1)]; % [sig f mu phi]
UppConstr = [space_sig(end) space_f(end) space_mu(end) space_phi(end)]; % [sig f mu phi]

Fs = 100;
dt = 1/Fs;
time = (-1:dt:1e2 + 10);
modelfun = @(sig, f, mu, phi) exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ); % original Kernel model
% costfun = @(p) sum( (Y' - fun(p(1),p(2),p(3))).^2 ); % Non Linear Least Square 
%% iterations
warning('off');
numIter = 100;
est = [];

for n = 1:numIter
        % draw 4 indx int from 1:200 and define parameters to estimate
        randIdx = randi(200,4,1);
        p = [space_sig(randIdx(1)) space_f(randIdx(4)) space_sig(randIdx(1))/2 + space_mu(randIdx(2)) space_phi(randIdx(3)) ];
        UppConstr(3) = space_sig(randIdx(1))/2 + space_mu(randIdx(2)); % [sig f mu phi]        
        fprintf('starting iteration number %d...\n', n);
        simulate4Parameters( p );
end

    function [] = simulate4Parameters( p )
        % Nested function for simulating the process for 4 new parameters
        
        % define LinearKernel
        sig_draw = p(1);       f_draw = p(2); 
        mu_draw = p(3);         phi_draw = p(4); % range abs [0.005, 2*Fs]
        targetfun = @(t) exp(-(t-mu_draw).^2/sig_draw^2).*sin( 2 * pi * f_draw *(t-phi_draw) );
        
%         figure(1);
%         plot(time,targetfun(time));
%         figure(2);plot(time,modelfun(p(1),p(2),p(3),p(4)))
%         figure(1);plot(time,Y_est)
%         pause(1);

        %% Optimization 
        % using the same optimization procedure as in section 2.1
        
        X = time;
        Y = targetfun(time);
        costfun = @(x) sum( (Y - modelfun(x(1),x(2),x(3),x(4))).^2 ); % SA with direct search  [sig f mu phi]
        
        % uses matlab's Simmulated annealing with pattern search
        % otimization
        x0 = [1 1 1 1]; % start point
        fprintf('starting optimization with SA ...\n');
        
        %% SA
        
        [p_est,~,~,output] = annealData2Acorr(x0, LowConstr,UppConstr,costfun);
        Y_est = modelfun( p_est(1), p_est(2), p_est(3), p_est(4)); % [sig f mu phi]
        numcoeff = 4;
        N = length(X);
        dfe = N - numcoeff; % degrees of freedom - number of observations - num of coefficients
        res = Y - Y_est;        
        %% Evaluating goodness of optimization and fit
        
        goodness = iGoodnessStructure(Y_est,[],res,dfe,N);
        optimResults = struct('X',X,'Data',Y,'estimation', Y_est,'optimOut', output);
        paramData = p; 
        paramEst = [ abs(p_est(1)), p_est(2), p_est(3), p_est(4) ];
        paramdiff = (paramData - paramEst)./mean([paramData;paramEst]);
        est_RMS = sqrt( (paramdiff * paramdiff' )/numel(paramdiff) ); % RMS of the estimation       
        PEresults = struct('GoodnessStatistics',goodness,'paramData',paramData,...
            'paramEst',paramEst,'estimationRMS', est_RMS);
        estIteration = struct('Results',optimResults,'PEresults',PEresults);
        est = [est estIteration];
    end

end


