%% Fake example for bias demonstration

% generate signal and white gaussian noise
x =  exp( -((-pi:0.01:pi) +1).^2/(1^2) ).*sin( 2*(-pi:0.01:pi) );
n = randn(1e6,1);

% find acorr of noise
maxlags = 300;
R_n = xcorr( n, maxlags,'unbiased');

% 'same' convolve with original signal and compare
x_conv = conv(x, R_n,'same' ); % x{k} is biased due to the 'same' property (the signal is not really WSS)

figure(1) 
plot(x); hold on;
plot(x_conv,'r'); hold off;
title('compare numerical bias with ''same'' conv')

figure(2) 
plot(conv( x, fliplr(x))); hold on;
plot(conv( x_conv, fliplr(x_conv)),'r'); hold off;
title('compare numerical bias with ''same'' conv')

% 'full' convolve with original signal and compare
x_conv = conv(R_n, x,'full'); % x{k} is biased due to the 'same' property (the signal is not really WSS)
x_conv = x_conv(1:length(x));

figure(2) 
plot(x); hold on;
plot(x_conv,'r'); hold off;
title('compare numerical bias with ''same')