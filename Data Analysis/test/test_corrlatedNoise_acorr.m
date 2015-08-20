%% test correlated noise structure

% create white gaussian noise and correlate
factor = 100;
w = 1*randn(1e6,1)';
w_correlated = repmat(w,[factor 1]);
w_correlated = w_correlated(:);

% find acorr
maxlags = 300;
[w_corr,lags] = xcorr(w_correlated, maxlags,'unbiased');

% create w_corr_tag
w_corr_tag = zeros(length(w_corr),1);
mainlobe = cumsum(w.*w);
mainlobe = mainlobe(end-factor:end);
mainlobe = [flipud(mainlobe) ; mainlobe(2:end)];
w_corr_tag = padarray(mainlobe, [(length(w_corr)-length(mainlobe))/2 0]);

% plot
figure();
plot(lags, w_corr,lags, w_corr_tag);
