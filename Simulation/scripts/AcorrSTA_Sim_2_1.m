function [est, proc] = AcorrSTA_Sim_2_1
%% set parameters and simualte process 
% define the numIter, the underlying WGN and the parameter space and
% simulate the process

% simulate process
Fs = 100;       N = 5e6;                        
noiseVar = 1;   noiseMean = 0;                          
process = PointProcessSynt( Fs, N, noiseMean, noiseVar );

% define parameter space
space_sig = linspace(1e-1, 1 , 200);       space_mu = linspace(0, 1e2, 200); 
space_fi = linspace(0, 2*pi, 200);         space_f = linspace(0, 8 ,200); % [rad] & [hz] range abs[0.005, 2*Fs]

% define upper and lower bounderies for optimization
constraintsLow = [eps 0 0]; % [ N f sig ]
constraintsUpp = [Inf space_f(end) space_sig(end)];

% define NlKernel and model for estimation
mu2 = 0;        sig2 = 0.8;
value1 = 'exp';                 value2 = @(t) exp( mu2 + sig2 * t);
NlKern = struct('type', value1, 'model', value2 ,'parameters', [mu2, sig2]);
process.NlinKernel = NlKern;
amp = 0.1;
modelEq = 'N * (exp((-x.^2)/(2*s^2)).*cos((2*pi*f).*x))-exp((-x.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0)';
fun = @(N,f,s) N * (exp((-X.^2)/(2*s^2)).*cos((2*pi*f).*X))-exp((-X.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0);
costfun = @(p) sum( (Y' - fun(p(1),p(2),p(3))).^2 ); % Non Linear Least Square 
    
%% iterations
warning('off');
numIter = 100;
[est, proc] = deal( cell(numIter,1) );

for n = 1:numIter
        
        randIdx = randi(200,4,1);
        p = [space_sig(randIdx(1)) space_mu(randIdx(2)) space_fi(randIdx(3)) space_f(randIdx(4))];
        
        % validates kernels position is positive   
        while (p(2) - (2 * sqrt(2*log(10)) * p(1))/2) <= 0
            randIdx = randi(200,4,1);
            p = [space_sig(randIdx(1)) space_mu(randIdx(2)) space_fi(randIdx(3)) space_f(randIdx(4))];
        end   
        
        fprintf('starting iteration number %d...\n', n);
        simulate4Parameters( p, n );
end

    function [] = simulate4Parameters( p, Iter )
        % Nested function for simulating the process for 4 new parameters
        
        % define LinearKernel
        sig1 = p(1);       mu1 = p(2); 
        fi = p(3);         f = p(4); % range abs [0.005, 2*Fs]
        value1 = 'gaussSine';           value2 = @(t) exp(-(t-mu1).^2/sig1^2).*sin( 2 * pi * f *(t-fi) );
        maxlags = ceil((2 * sqrt(2*log(10)) * sig1)/process.dt); % the Full Width a 0.1 of Max  
        support = ( ceil(mu1/process.dt - maxlags/2):ceil(mu1/ process.dt + maxlags/2) )';
        timelag = (-maxlags : maxlags) .* process.dt;
        LinKern = struct('type', value1, 'model', value2, 'support', support,...
            'maxlags',maxlags,'timelag', timelag, 'parameters', [sig1, mu1, pi, f]);
        process.linKernel = LinKern;
        
        % define FIR Window lowpass filter for a preproceesing step of the estimation
        fc = space_f(end) + 5;
        Wn = (2/Fs)*fc; 
        if 2*maxlags <= 300 % data must be more than 3 times the filter order
            ord = floor(maxlags*2/3)-1;
        else
                ord = 100;
        end
        filt_coeff = fir1(ord,Wn,'low',kaiser(ord+1,1));
        %% Generate process and estimate parameters from simulated process

        process = process.generateProcess( amp );
        estimation = PointProcessSynt( Fs, N, noiseMean, noiseVar );
        estimation.spiketimes = process.spiketimes;
        estimation = estimation.estimatesStatsProperties( process.linKernel.maxlags, process.acorr.lambdaCorr );
        timeAxis = (-process.linKernel.maxlags : process.linKernel.maxlags) .* process.dt;
        estimation.acorr.lags = timeAxis;
        %% Parameter estimation
        % Estimation of 'f' and 'sig1' parameters out of the correlation function

        est{Iter} = estimation.paramEstim( process, timeAxis , constraintsLow, constraintsUpp, filt_coeff);
        proc{Iter} = process;
    end

end

