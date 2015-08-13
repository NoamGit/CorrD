%% Negative auto-correlations for counting process

% Load data and take spike samples
D = load('SpikeTimeGauss1B.mat'); % change path
sr = 1e4; % sampling rate 10 [kHz]
dt = 1/sr;
dataIdx = 5; % choose index between [1-23]
spike_samples = D.TT(dataIdx).sp;
spike_times = spike_samples .* dt;

% create counting process and time axis
Tmax = max(spike_times);
Tmin = min(spike_times);
time_axis = (Tmin:dt:Tmax);
counting_process = histc(spike_times, time_axis);

% compute correlations and plot
maxlags = 1e3;
bias = 'unbiased';
[R, lag] = xcorr(counting_process, maxlags, bias);
R = R./dt^2;

% plot
figure();
plot( lag.*dt ,R );
title('acorr example - single data instance');
assert( ~any(R < 0 ) , 'Attention: Not all values are positive');