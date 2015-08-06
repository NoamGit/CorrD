clear all; close all; clc
% TestRetrievalConcept
%
% Description: Obtain L filter parameters out of syntesized data
% This file demonstrates the backprocedure of estimating the parameters of
% the linear filter out of the AC of syntesized data.
% Result: works well for a certain range of parameters
%% Syntesize Counting Processes

GWN = 0 + 3 * randn(1e6,1); % Gaussian white noise
parameters = [0.5 5 -1.2 0]'; % parameters for kernel - sigma f mu phi
dt = 1e-2; % dt = 1/N
range_h = (-3:dt:-0); % range of kernel (we need something short)
h = @(x) exp(-(range_h-x(3)).^2/x(1)^2).*sin(x(2)*(range_h-x(4))); % kernel
CGP = conv( GWN, 0.1*h(parameters),'valid')'; % Correlated Gaussian process
Int = exp(CGP)'; % Intensity function exponential
t = (0:length(Int)-1)'*dt; % time
spiketimes = simpp(Int, dt);
maxlags = 500; % specifing max lags for efficency
Int_AC = xcorr(Int, maxlags, 'unbiased'); % intensity AC 
valueTau0 = Int_AC(ceil(length(Int_AC)/2)); % value at tau == 0
fprintf(['E(lambda) = ', num2str(length(spiketimes)/t(end)),' Hz \n']); % print mean firing rate for checking feasibility    
%% Back procedure - finding freq and sigma (parameters) out of the AC of PP

Int_AC_retrieved = ppcorr(spiketimes, dt, maxlags, 'unbiased', 1); %CP calculation in new resolution using michaels function ppcorr
E_lambda = ppmean(spiketimes); % mean firing rate calculation using MK function
tau = dt*(-maxlags:maxlags)';
tau0_indx =  ceil(length(Int_AC)/2);
Int_AC_retrieved(tau0_indx) = valueTau0; % assuming that we know tau = 0
set(gca,'FontSize',14);
plot([Int_AC_retrieved Int_AC]); % Intensity AC comparsion
title('Intesity back and forward procedure');
legend('retrieved', 'original'); axis tight;
%% Comparing AC of STA and AC of retrieved CGP for NL = exp

sigma = sqrt(log(Int_AC_retrieved(tau0_indx)/E_lambda.^2)); % estimated sigma exp 
CGP_AC = (1/sigma.^2)*log(Int_AC_retrieved/E_lambda.^2);% exp 
STA_AC = xcorr(conv(h(parameters),GWN), maxlags,'unbiased'); % CGP AC is actually equivalent to STA AC (stimulus is GWN)

% PLOTS
figure();set(gca,'FontSize',14); 
% plot(tau, (CGP_AC)./std(CGP_AC),tau ,(STA_AC)./std(STA_AC))
plot(tau, normax(CGP_AC),tau ,normax(STA_AC))
title(['AC comparsion of CGP and STA ', '\sigma = ' num2str(parameters(1)), ' \freq = ' num2str(parameters(2))]);
legend('norm GCP AC','norm STA AC')
xlabel('time [sec]'); ylabel('norm Amp');
%% Parameter estimation using the developed analytical AC

Weights = 0.2 * ones(length(tau),1); % weight function for interpolation 
Weights(tau0_indx-80:tau0_indx+80) = 0.8; Weights(tau0_indx) = 0.1;
CGP_AC_sm = feval(createSmoothing(CGP_AC, Weights),(1: size(CGP_AC,1)));
weigths = calcWeights(CGP_AC_sm);
[CGP__AC_fit,gof] = createFit2(tau, CGP_AC, weigths);
[xData, yData, weights] = prepareCurveData( tau, CGP_AC,weigths);
% Plot fit with data.
handl = plot( CGP__AC_fit, xData, yData );
legend( handl, 'GCP data', 'Parametric estimation', 'Location', 'NorthEast' );
%% 

STA = compute_sta(CGP', spiketimes, length(h(parameters)), dt);

figure();
STA_AC_2 = xcov( STA, maxlags,'unbiased')'; % STAs maxlags array size
plot(tau, (STA_AC_2)./std(STA_AC_2),tau ,(STA_AC)./std(STA_AC))
title('AC comparsion of CGP_AC and STA_AC ');
legend('norm GCP AC','norm STA AC')
xlabel('time [sec]'); ylabel('norm Amp');
