clear all; close all; clc
% TestRetrievalConcept2
%
% Description: This is similar to TestRetrievalConcept1, but is inspired from 
% AcorrV13. Obtain L filter parameters out of syntesized data
% This file demonstrates the backprocedure of estimating the parameters of
% the linear filter out of the AC of syntesized data.
%% Defining system's variabels

Fs = 100;                         %100Hz sample rate
dt = 1/Fs;
N = 1e6;                            %num of sampels
T = N/Fs;                           %signal length
t = linspace(0,T,N);                %timeline % print mean firing rate for checking feasibility 
sigmaP = 5;                             %arbitrary sigma [0.5,3]
f = 0.6;                                 %arbitrary f limited to abs[0.005,2Fs]
phi = 4;                                %arbitrary fi
mu = 70;                               %arbitrary mu [8,13]
maxlags = ceil((2*sqrt(2*log(10))*sigmaP)/dt); %the Full Width a tenth of maximum (FW0.1M), this represent the relevat data when we trim our signals               CONSTRAIN!!                                      
% maxlags = 500;
% range_h = (mu/dt-maxlags/2) : (mu/dt + maxlags/2);
% h =@(x) exp(-(range_h-x(3)).^2/x(1)^2).*sin(x(2)*(range_h-x(4))); % kernel
% h_trim = h([sigma f mu phi]);
% 
h = exp(-(t-mu).^2/sigmaP^2).*sin(f*(t-phi));
h_trim = h((mu/dt-maxlags/2) : (mu/dt + maxlags/2));
%System's Input - white gaussian noise (WGN) mean=0, var=1
GWN = 0 + randn(N, 1);
CGP_AC_F = xcorr(conv(h_trim,GWN), maxlags,'unbiased'); % CGP AC
%% Syntesize Counting Processes

Amp = 0.1;                      %CONSTRAIN!!!
CGP = conv(GWN,Amp * h_trim,'same');
sigmaI = 0.8; muI = 0;
Int = exp( muI + sigmaI*(CGP) ); 
spikestimes = simpp(Int,dt); 
%% Back procedure - finding freq and sigma (parameters) out of the AC of PP

Int_AC =  ppcorr(spikestimes, dt, maxlags, 'unbiased', 1);
[R_I,tau2] = xcorr(Int, maxlags,'unbiased');
tau_0 = find(tau2 == 0);
Int_AC_retrieved = Int_AC;
Int_AC_retrieved(tau_0) = R_I(tau_0); % assuming that we know tau = 0
idx_tau_0 =  ceil(length(Int_AC)/2);         
E_lambda = ppmean(spikestimes);

sigma = sqrt(log(Int_AC_retrieved(idx_tau_0)/(E_lambda.^2)));
CGP_AC_B = (1/sigma^2)*log(Int_AC_retrieved/(E_lambda.^2));             %according to michael cacls in p.20

% PLOTS
figure();set(gca,'FontSize',14); 
plot( tau2, CGP_AC_B/norm(CGP_AC_B),'-g', tau2, CGP_AC_F/norm(CGP_AC_F))
% title(['AC comparsion of CGP and STA ', '\sigma = ' num2str(parameters(1)), ' \freq = ' num2str(parameters(2))]);
legend('norm GCP AC','norm STA AC')
xlabel('time [sec]'); ylabel('norm Amp');
axis tight;