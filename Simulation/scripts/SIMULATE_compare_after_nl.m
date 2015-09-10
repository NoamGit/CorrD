% Complete simulation of all important functions and comparison of acorr
% after the non linearity

% Part 1 - Generate Syntetic Data
%% generate stimulus
dt_stim = 0.01;
N = 1e6;
noiseVar = 1;
noiseMean = 0;
stim = noiseMean + noiseVar * randn(N,1);
t_stim = linspace(0,N*dt_stim,N);

c = 10;
syntData = [];
%% create all processes
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
    amplitude = 0.1;
    model_kernel = @(t) amplitude * exp(-(t-p_kernel(2)).^2/p_kernel(1)^2).*sin( p_kernel(4) *(t-p_kernel(3)) );
    t_kernel = linspace(-30*dt_kernel, 0*dt_kernel, 1/dt_stim); 
    
%     figure(1);plot(t_kernel, model_kernel(t_kernel));% show kernel
    %% create nl kernel
    
    theta1_spc = linspace(0.0001, 0.0004, 200);       
    theta2_spc = linspace(0.1, 0.3, 200); 
    theta3_spc = linspace(2,5, 200);
    randIdx = randi(200,4,1);
    % draw values and find kernel
    p_nl = [theta1_spc(randIdx(1)) theta2_spc(randIdx(2)) theta3_spc(randIdx(3))];
    model_nl = @(x) p_nl(1)* x.^2 + p_nl(2) * x + p_nl(3);  % define parameter space for nl
%     model_nl = @(x) exp(x);
    %% create CGP, lambda and spiketimes
    
    amp = 1e3;
    CGP = conv( stim,amp * model_kernel(fliplr(t_kernel)),'full');
    lambda = model_nl(CGP); 
    spiketimes = simpp( lambda , dt_stim);
    
    nl_domain = (min(CGP):dt_stim:max(CGP));
    nl_image = model_nl(nl_domain);

%     figure(2); plot(nl_domain, nl_image);% show nl
%     CountigProcess = histc(spiketimes,t_stim);
%     R_dN = xcov(CountigProcess  , maxlags, 'unbiased' )/(dt_stim^2);
%     [R_lambd,tau] = xcov(lambda, maxlags,'unbiased');
%     R_dN(maxlags) = [];R_lambd(maxlags) = [];tau(maxlags) = [];
%     plot(tau*dt_stim, normax(R_dN)); hold on;
%     plot(tau*dt_stim,normax(R_lambd)); hold off;
%     %% compare sta
%     
%     sta = compute_sta( stim, spiketimes, length(t_kernel), 0 , dt_stim);
%     plot(t_kernel,normalize(sta,6),t_kernel,normalize(model_kernel(t_kernel),6));
    %% store data
    
    kernel = struct('time',t_kernel,'val', model_kernel(fliplr(t_kernel)),'model',model_kernel,'param',p_kernel);
    nonLinkernel = struct('domain',nl_domain,'image',nl_image,'model',model_nl,'param',p_nl);
    datainstance = struct('linKernel',kernel,'nonlinKernel',nonLinkernel,...
        'CGP',CGP,'lambd',lambda,'spiketimes',spiketimes);  
    syntData = [syntData  datainstance];
    display(['finished syntesizing process ',num2str(n)])
end
    %% store generated data in process object 
    
    % Load data to PointProcessData structure
    process = LoadDataIntoProcess( 'synt', dt_stim, syntData, stim );
    kernNAME = 'gaussSine';
    maxlags = ceil((10 * sqrt(2*log(10)) * max(space_sig))/process.dt); 
    process.maxlags = maxlags;
    % Part 2 - test custom made functions
    %% try calc STA
    
    linKernel_orig = cellfun(@(x) x.val, {syntData.linKernel},'UniformOutput', false);
    process = process.CalcSTA(length(linKernel_orig{1}),0, 'original stimulus'); 
    % plot
    plotCompare({process.STA.realSTA},linKernel_orig, t_kernel, 10 ,'STA of cell ',6,(1:10))
    legend('STA','Kernel');
    %% try nlestimation
    
    process.CGP = cellfun(@(x) conv( stim,amp * x.model(fliplr(t_kernel)),'same'),{syntData.linKernel},'UniformOutput',false);
    process.CP = arrayfun(@(x) histc(process.spiketimes{x}, t_stim),(1:c),'UniformOutput',false);
    poly_order = 2;
    process = process.nlestimation(poly_order, 0, amp, 'resampled stimulus','Synt Data'); % stimulus has same dt as procees -> resampled
    % plot
    nl_domain = cellfun(@(x) x.domain, {syntData.nonlinKernel},'UniformOutput', false);
    nl_domain_norm = cellfun(@(x) x/std(x), nl_domain,'UniformOutput',false);
    nl_image = cellfun(@(y,x) y.model(x), {syntData.nonlinKernel},nl_domain,'UniformOutput', false);
    fullmodel = @(theta,x)  theta(1)* x.^2 + theta(2) * x + theta(3);
    nl_est_image = cellfun(@(c,x) fullmodel(c.polyfit,x), {process.NlinKernel.estimation},nl_domain,'UniformOutput', false);
    plotCompare(nl_image,nl_est_image,nl_domain, 10 ,'nl of cell ', 6 ,(1:10));
    legend('est NL','NL');
    %% try compareInNlPlane


