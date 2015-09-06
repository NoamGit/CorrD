% Complete simulation of all important functions and comparison of acorr
% after the non linearity

% Part 1 - Generate Syntetic Data
%% generate stimulus
dt_stim = 0.01;
N = 1e6;
noiseVar = 1;
noiseMean = 0;
stim = noiseMean + noiseVar * randn(N,1);
t_stim = (0:dt_stim:N*dt_stim);

c = 10;
syntData = [];
for n = 1:c % build c processes
    %% create linear kernels

    % define parameter space for linear kernel
    space_sig = linspace(0.1, 0.15, 200);       space_mu = linspace(-0.4, -0.1, 200); 
    space_fi = linspace(pi, 2*pi, 200);         space_f = linspace(8, 15 ,200); % [rad] & [hz] range abs[0.005, 2*Fs]
    randIdx = randi(200,4,1);
    % draw values and find kernel
    p_kernel = [space_sig(randIdx(1)) space_mu(randIdx(2)) space_fi(randIdx(3)) space_f(randIdx(4))];
    dt_kernel = 1/30;
    amplitude = 10;
    model_kernel = @(t) amplitude * exp(-(t-p_kernel(2)).^2/p_kernel(1)^2).*sin( p_kernel(4) *(t-p_kernel(3)) );
    t_kernel = linspace(-30*dt_kernel, 10*dt_kernel, 1/dt_stim); 
    
%     figure(1);plot(t_kernel, model_kernel(t_kernel));% show kernel
    %% create nl kernel
    
    % define parameter space for nl
    theta1_spc = linspace(0.001, 0.004, 200);       
    theta2_spc = linspace(0.0005, 0.001, 200); 
    theta3_spc = linspace(2,5, 200);
    randIdx = randi(200,4,1);
    % draw values and find kernel
    p_nl = [theta1_spc(randIdx(1)) theta2_spc(randIdx(2)) theta3_spc(randIdx(3))];
    model_nl = @(x) p_nl(1)* x.^2 + p_nl(2) * x + p_nl(3);
    %% create CGP, lambda and spiketimes
    
    amp = 0.01;
    CGP = conv(stim ,amp * model_kernel(t_kernel),'same');
    lambda = model_nl(CGP); 
    spiketimes = simpp( lambda , dt_stim);
    
    nl_domain = (min(CGP):dt_stim:max(CGP));
    nl_image = model_nl(nl_domain);
    
%     figure(2); plot(nl_domain, nl_image);% show nl   
    %% store data
    
    kernel = struct('time',t_kernel,'val', model_kernel(t_kernel),'model',model_kernel,'param',p_kernel);
    nonLinkernel = struct('domain',nl_domain,'image',nl_image,'model',model_nl,'param',p_nl);
    datainstance = struct('linKernel',kernel,'nonlinKernel',nonLinkernel,...
        'CGP',CGP,'lambd',lambda,'spiketimes',spiketimes);  
    syntData = [syntData  datainstance];
end
    %% store generated data in process object 
    
    % Load data to PointProcessData structure
    process = LoadDataIntoProcess( 'synt', dt_stim, syntData, stim );
    kernNAME = 'gaussSine';
    maxlags = ceil((0.621 * sqrt(2*log(10)) * 1)/process.dt); 
    process.maxlags = maxlags;
    % Part 2 - test custom made functions
    %% try calc STA
    
    linKernel_orig = cellfun(@(x) x.val, {syntData.linKernel},'UniformOutput', false);
    process = process.CalcSTA(length(linKernel_orig{1}),10, 'original stimulus'); 
    plotCompare({process.STA.realSTA},linKernel_orig, t_kernel, 4 ,'STA of cell ',0,(1:10))
    %% try nlestimation
    
    poly_order = 2;
    process = process.nlestimation(poly_order, showFlag, 'original stimulus'); 
    %% try compareInNlPlane


