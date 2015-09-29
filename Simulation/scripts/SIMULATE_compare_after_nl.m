% Complete simulation of all important functions and comparison of acorr
% after the non linearity

% Part 1 - Generate Syntetic Data
%% generate stimulus
dt_stim = 0.01;
N = 1e5;
noiseVar = 1;
noiseMean = 0;
stim = noiseMean + noiseVar * randn(N,1);
t_stim = linspace(0,N*dt_stim,N);

c = 6;
nonLinMethod = 'exp';
% nonLinMethod = 'poly';
%% create all processes
syntData = [];
for n = 1:c % build c processes
    %% create linear kernels

    % define parameter space for linear kernel
    space_sig = linspace(0.1, 0.15, 200);       space_mu = linspace(-0.6, -0.2, 200); 
    space_fi = linspace(pi, 2*pi, 200);         space_f = linspace(8, 15 ,200); % [rad] & [hz] range abs[0.005, 2*Fs]
    randIdx = randi(200,4,1);
    % draw values and find kernel
    p_kernel = [space_sig(randIdx(1)) space_mu(randIdx(2)) space_fi(randIdx(3)) space_f(randIdx(4))];
    maxlags = 2*ceil((2*sqrt(2*log(10))*space_sig(randIdx(1)))/dt_stim);   
    dt_kernel = 1/30;
    amp_kern = 1;
    model_kernel = @(t) amp_kern * ( exp(-(t-p_kernel(2)).^2/p_kernel(1)^2).*sin( p_kernel(4) *(t-p_kernel(3)) ) );
    t_kernel = linspace(-30*dt_kernel, 0*dt_kernel, 1/dt_stim);   
    %% create nl kernel
    theta1_spc = linspace(0.00005, 0.00025, 200);       
    theta2_spc = linspace(0.05, 0.3, 200); 
    theta3_spc = linspace(170,200, 200);
    theta4_spc = linspace(18, 20, 200);       
    theta5_spc = linspace(0.165,0.20, 200); 
        
    h = tic;
    readyFlag = 1;
    while(readyFlag && toc(h) < 60)
        randIdx = randi(200,5,1);
        % draw values and find kernel
        p_nl = [theta1_spc(randIdx(1)) theta2_spc(randIdx(2)) theta3_spc(randIdx(3)) theta4_spc(randIdx(4)) theta5_spc(randIdx(5))];
        switch nonLinMethod
            case 'poly'
                amp_stim = 800;
                model_nl = @(x) p_nl(1)* x.^2 + p_nl(2) * x + p_nl(3);  % define parameter space for nl
            case 'exp'
                amp_stim = 5.58;
                model_nl = @(x) p_nl(4) * exp(x.* p_nl(5));  
        end
        %% create CGP, lambda and spiketimes

        CGP = conv( stim,amp_stim * model_kernel(fliplr(t_kernel)),'full');
        CGP = CGP(1:end-length(t_kernel)+1);
        lambda = model_nl(CGP); 
        spiketimes = simpp( lambda , dt_stim);

        nl_domain = (min(CGP):dt_stim:max(CGP));
        nl_image = model_nl(nl_domain);
        if ~any( nl_image <0 ) &&...
                numel(spiketimes)/N > 5 &&...
                numel(spiketimes)/N < 20
            readyFlag = 0;
        end
    end
    %% store data
    
    kernel = struct('time',t_kernel,'val', model_kernel(t_kernel),'model',model_kernel,'param',p_kernel);
    nonLinkernel = struct('domain',nl_domain,'image',nl_image,'model',model_nl,'param',p_nl);
    datainstance = struct('linKernel',kernel,'nonlinKernel',nonLinkernel,...
        'CGP',CGP,'lambd',lambda,'spiketimes',spiketimes);  
    syntData = [ syntData  datainstance ];
    display(['finished syntesizing process ',num2str(n)])  
    %% Validate data
    
    assert(~any( nl_image < 0 ), 'Non linearity is not def positiv');
    assert(numel(spiketimes)/N > 5 && numel(spiketimes)/N < 20, 'mean firing rate is to low or too high');
    %% iteration plot
    
%     figure(1);plot(t_kernel, amp_stim * model_kernel(t_kernel));% show kernel
% 
%     figure(2); plot(nl_domain, nl_image);% show nl
%     
%     CountigProcess = histc(spiketimes,t_stim);
%     R_dN = xcov(CountigProcess  , maxlags, 'unbiased' )/(dt_stim^2);
%     [R_lambd,tau] = xcov(lambda, maxlags,'unbiased');
%     R_dN(maxlags) = [];R_lambd(maxlags) = [];tau(maxlags) = [];    
%     figure(3);
%     plot(tau*dt_stim, normax(R_dN)); hold on;
%     plot(tau*dt_stim,normax(R_lambd),'r'); hold off;
%     
%     figure(4)
%     sta = compute_sta( stim, spiketimes, length(t_kernel), 0 , dt_stim);
%     plot(t_kernel,normalize(sta,6),t_kernel,normalize(model_kernel(t_kernel),6));
% 
%     numel(spiketimes)/N
end
    %% store generated data in process object 
    
    % Load data to PointProcessData structure
    process = LoadDataIntoProcess( 'synt', dt_stim, syntData, stim );
    maxlags = ceil((10 * sqrt(2*log(10)) * max(space_sig))/process.dt); 
    process.maxlags = maxlags;
    % Part 2 - test custom made functions
    %% check calc STA
    
    linKernel_orig = cellfun(@(x) amp_stim * x.val, {syntData.linKernel},'UniformOutput', false);
    process = process.CalcSTA(length(linKernel_orig{1}), 0, amp_stim, 'original stimulus'); 
    % normalize STA 
    sta_norm = cellfun(@(x, y) normalize(x,6).*std(y)+mean(y) ,{process.STA.realSTA},linKernel_orig,'UniformOutput',false);
    for k = 1:c
        process.STA(k).realSTA = sta_norm{k};
    end
    
    % plot
    plotCompare({process.STA.realSTA},linKernel_orig, t_kernel, c ,'STA of cell ',0,(1:c))
    legend('STA','Kernel');
    %% check nlestimation
    
    cgp_cellArray = cellfun(@(x) conv( stim,amp_stim * x.model(fliplr(t_kernel)),'full'),{syntData.linKernel},'UniformOutput',false);
    process.CGP = cellfun(@(x) x(1:end-length(t_kernel)+1) ,cgp_cellArray,'UniformOutput',false);
    process.CP = arrayfun(@(x) histc(process.spiketimes{x}, t_stim),(1:c),'UniformOutput',false);
    
    switch nonLinMethod
        case 'poly'
            poly_order = 2;
            process = process.nlestimation(poly_order, 0, amp_stim,'poly', 'resampled stimulus','Synt Data'); % stimulus has same dt as procees -> resampled
            nl_domain = cellfun(@(x) downsample(x.domain,5e4), {syntData.nonlinKernel},'UniformOutput', false);
            nl_image = cellfun(@(y,x) y.model(x), {syntData.nonlinKernel},nl_domain,'UniformOutput', false);
            fullmodel = @(theta,x)  theta(1)* x.^2 + theta(2) * x + theta(3);
        case 'exp'
            poly_order = 0;
            process = process.nlestimation(poly_order,0, amp_stim,'exp', 'resampled stimulus','Synt Data'); % stimulus has same dt as procees -> resample
            nl_domain = cellfun(@(x) downsample(x.domain,10), {syntData.nonlinKernel},'UniformOutput', false);
            nl_image = cellfun(@(y,x) y.model(x), {syntData.nonlinKernel},nl_domain,'UniformOutput', false);
            fullmodel = @(theta, x) theta(1) * exp(x.* theta(2));
    end
    % plot
            nl_est_image = cellfun(@(c,x) fullmodel(c.fitParam,x), {process.NlinKernel.estimation},nl_domain,'UniformOutput', false);
            plotCompare(nl_image,nl_est_image,nl_domain, c ,'nl of cell ', 6 ,(1:c));
            legend('est NL','NL');   
    %% try comparing RdN 
    
    flagDoublet = 0; flagCompare = 0;
    process = process.PreProcessCP( 'estimatedNL', flagDoublet, flagCompare );
    compareInNlPlane( process );