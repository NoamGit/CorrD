rng default

% create sine
Fs = 1000;
t = linspace(0,1,Fs);
N = length(t);
x = cos(2*pi*100*t);
y = x + 1.5*randn(size(t));

% estimate PSD
Y = fft(y);
Pyy = Y.*conj(Y)/N;
f = Fs/N*(0:N/2);
figure;
plot(f,Pyy(1:N/2+1))
title('Power spectral density')
xlabel('Frequency (Hz)')

% FIR Window lowpass filter designed using the FIR1 function.
fc = 120;
Wn = (2/Fs)*fc;
b = fir1(90,Wn,'low',kaiser(91,1));
fvtool(b,1,'Fs',Fs)

%apply filter an visualize
y1 = filtfilt(b,1,y);
plot(t(1:100),x(1:100))
hold on
plot(t(1:100),y1(1:100),'g')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original Signal','Filtered Data')

% estimate PSD of filtered signal
Y1 = fft(y1);
Pyy1 = Y1.*conj(Y1)/N;
figure;
plot(f,Pyy1(1:N/2+1))
title('Power spectral density')
xlabel('Frequency (Hz)')