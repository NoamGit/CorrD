% main
% clear all; close all; clc
%% set parameters and simualte process 

% simulate process
Fs = 100;       N = 1e6;                        
noiseVar = 1;   noiseMean = 0;                          
process = PointProcessSynt( Fs, N, noiseMean, noiseVar );

% define LinearKernel
sig1 = 3;        mu1 = 70; % [sec]
pi = 4;         f = 3; % [rad] & [hz] range abs[0.005, 2*Fs]

% validates kernels positive position 
% if mu1 - (2 * sqrt(2*log(10)) * sig1)/2 <= 0
%     continue;
% end
value1 = 'gaussSine';           value2 = @(t) exp(-(t-mu1).^2/sig1^2).*sin( f *(t-pi) );
maxlags = ceil((2 * sqrt(2*log(10)) * sig1)/process.dt); % the Full Width a 0.1 of Max in [samp]
support = ( ceil(mu1/process.dt - maxlags/2):ceil(mu1/ process.dt + maxlags/2) )';
timelag = (-maxlags : maxlags) .* process.dt;
LinKern = struct('type', value1, 'model', value2, 'support', support,...
    'maxlags',maxlags,'timelag', timelag, 'parameters', [sig1, mu1, pi, f]);

% define NlKernel set kernels 
mu2 = 0;        sig2 = 0.8;
value1 = 'exp';                 value2 = @(t) exp( mu2 + sig2 * t);
NlKern = struct('type', value1, 'model', value2 ,'parameters', [mu2, sig2]);
process.linKernel = LinKern;
process.NlinKernel = NlKern;
amp = 0.1;
%% Generate process and estimate parameters from simulated process

process = process.generateProcess( amp );
estimation = PointProcessSynt( Fs, N, noiseMean, noiseVar );
estimation.spiketimes = process.spiketimes;
estimation = estimation.estimatesStatsProperties( process.linKernel.maxlags, process.acorr.lambdaCorr );
%% Plotting  
% Graphical validation of the estimated acorr versus the acoor of the generated GCP
% empirically, normalizing by 3 yields good resutls

figure();set(gca,'FontSize',14); 
plot( timelag, estimation.acorr.CGP_Corr./norm(estimation.acorr.CGP_Corr,3) ); hold on;
plot( timelag ,process.acorr.CGP_Corr./norm(process.acorr.CGP_Corr,3),'-g' ); hold off;
legend('normalized GCP acorr process','normalized GCP acorr reconstruction')
xlabel('time [sec]'); ylabel('norm Amp');
axis tight;
%% Parameter estimation
% Estimation of 'f' and 'sig1' parameters out of the correlation function

constraintsLow = [-Inf 0 0]; % [ N f sig ]
constraintsUpp = [Inf 3 50];
modelEq = 'N * (exp((-x.^2)/(2*s^2)).*cos(f.*x))-exp((-x.^2-f^2*s^4)/(2*s^2))';
Y = normax( estimation.acorr.CGP_Corr );
X = timelag;
estimation = estimation.paramEstim( process, timelag, modelEq, constraintsLow, constraintsUpp );

estimation.PEresults
estimation.PEresults.fitResults
estimation.linKernel.timelag = process.linKernel.timelag;
estimation.visualizePE( estimation )
estimation.visualizePE( process )
