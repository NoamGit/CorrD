%% Test SPMD

parallProcesses_NUM = 20;
numIterD = distributed(1:parallProcesses_NUM);
X = magic(10);
fun = arrayfun(@(n) (@(x) sum(n * x(:)')),(1:parallProcesses_NUM),'UniformOutput',false); 
out = zeros(parallProcesses_NUM,1);

spmd
        numIterLP = getLocalPart ( numIterD );
        for n = numIterLP
            
        fprintf('starting iteration number %d...\n', n);    
        % SA
        % CONTINUE - SPMD DOESNT WORKS WITH ANYNIMUS FUNC
        fun_ITER = fun{n};
        out(n) = fun_ITER(X);
        end
end

Composite2array( out )
%% Single optimization check

LowConstr = [eps -100 -100 0]; % [sig f mu phi]
UppConstr = [10 100 100 pi]; % [sig f mu phi]
x0 = [1 1 1 1];
numvar = 4; % for ga
PopulSize = 200;

    %% bad indexes [1 4 9 16]
    indx = 4;
    modelfun = @(sig, f, mu, phi) normax( exp((-(time-mu).^2)/sig.^2) .* sin( (2*pi*f).* (time - phi) ) ); % with normax
    costfun_TEMP = (@(x)sum((normax(realSTA{indx})'-modelfun(x(1),x(2),x(3),x(4))).^2));
    [p_est_TEMP,~,~,output_TEMP] = ga_Data2STA(numvar, LowConstr,UppConstr,costfun_TEMP);
%     [p_est_TEMP,~,~,output_TEMP] = gps_Data2STA(x0 , LowConstr,UppConstr,costfun_TEMP);
%     costfun_TEMP(p_est_TEMP)
    estSTA_TEMP = modelfun( p_est_TEMP(1), p_est_TEMP(2), p_est_TEMP(3), p_est_TEMP(4)); % [sig f mu phi]   
    res_TEMP = normax(realSTA{indx}) - estSTA_TEMP';  
    goodness_TEMP = iGoodnessStructure(estSTA_TEMP,[],res_TEMP,dfe,N);
%     plotFit( time, normax(realSTA{indx}) , estSTA_TEMP, 'STA 1 fit')
    disp(goodness_TEMP);
    
    %% using grid search
    % RESULTS - inefficient for 4 parameters retrival.
    
    NumDiv =50; % array of number of divisions along each dimension
    MinDeltaX = 1e-2; % array of minimum search values for each variable
    Eps_Fx = 1e-5; % tolerance for difference in successive function values
    MaxIter = 1e3; % maximum number of iterations
    indx = 19;    
%     costfun_grid = (@(x1,x2,x3,x4)sum((normax(realSTA{indx})'-modelfun(x1,x2,x3,x4)).^2));
    costfun_grid = @(x)sum((normax(realSTA{indx})'-modelfun(x(1),x(2),x(3),x(4))).^2);

    [XBest,BestF,Iters] = Grid_Search(numvar, LowConstr, UppConstr, NumDiv, MinDeltaX, Eps_Fx, MaxIter, costfun_grid);
    
    