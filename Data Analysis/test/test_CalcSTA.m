%% Test STA calculation

% load Segev's data
D = LoadDataIntoProcess( 1 ); 
DataArray = D.Data;

% turning the stimulus to white noise
stim = DataArray.stimulus.Yun;
% sTHRESCORR = 5.346e4; % stimulus threshold correlation
% stim = DataArray.stimulus.Yraw(1:sTHRESCORR);
% stim = stim - mean(stim);
% DataArray.stimulus.Yraw = stim;

% check spectrogram and PSD
figure
subplot(221)
spectrogram(stim,512,100)
title('Stimulu''s Spectrogram');
subplot(222)
plot(xcorr(stim,500,'unbiased'));
axis tight
title('Stimulus aCORR');
subplot(223)
pwelch(stim);
title('Power spectral density By Welch method')

% find STA
DataArray = DataArray.CalcSTA(35, 5); % 35 window size ; 5 left window edge in the STA
sta_SEGEV = RonenSTA(DataArray.spiketimes{6} * 1/DataArray.dt, stim, 1/DataArray.dt);

% plot results
plot(DataArray.STA.realSTA{3})
