% Example for PSD estimation for the process and estimation created in the
% NLP sheme

Fs = 1000;
t = linspace(0,1,Fs);
N = length(t);
x = exp(-(t-(300*1/Fs)).^2/(120*1/Fs)^2) .* cos(2*pi*10*t);
y = conv(x,5*randn(size(t)),'same');
Ry = xcorr(y,'unbiased');
Rx = xcorr(x,'unbiased');

% estimate PSD

Y = fft(Ry);
Pyy = Y.*conj(Y)/N;
f = Fs/N*(0:N/2);
figure;
plot(f,Pyy(1:N/2+1))
title('Power spectral density')
xlabel('Frequency (Hz)')
%% applying on simulated data proc
% you need to creat some data first with [est, proc] = AcorrSTA_Sim_2_1

ii =  4;
f_data = proc{ii}.linKernel.parameters(4);    
t = proc{ii}.t;
Fs = proc{ii}.Fs; 
N = 2 * proc{ii}.linKernel.maxlags;
x = proc{ii}.CGP;
Ry = proc{ii}.acorr.CGP_Corr;

% estimate PSD
Y = fft(Ry);
Pyy = Y.*conj(Y)/N;
f = Fs/N*(0:N/2);
figure;
plot(f,Pyy(1:N/2+1))
title(['Power spectral density. main freq - ',num2str(f_data)])
xlabel('Frequency (Hz)')
%% applying on simulated data est
% you need to creat some data first with [est, proc] = AcorrSTA_Sim_2_1

ii =  2;
f_data = proc{ii}.linKernel.parameters(4);    
t = est{ii}.t;
Fs = est{ii}.Fs; 
N = 2 * proc{ii}.linKernel.maxlags;
x = est{ii}.CGP;
Ry = est{ii}.acorr.CGP_Corr;
Rx = proc{ii}.acorr.CGP_Corr;

% estimate PSD
Y = fft(Ry);
Pyy = Y.*conj(Y)/N;
f = Fs/N*(0:N/2);
figure;
plot(f,Pyy(1:N/2+1))
title(['Power spectral density. main freq - ',num2str(f_data)])
xlabel('Frequency (Hz)')
%% filter est

% FIR Window lowpass filter designed using the FIR1 function.
fc = 30+2;
Wn = (2/Fs)*fc;
b = fir1(100,Wn,'low',kaiser(101,1));
% fvtool(Rx,1,'Fs',Fs)

%apply filter and visualize
Ry1 = filtfilt(b,1,Ry);

Y = fft(Ry1);
Pyy1 = Y.*conj(Y)/N;
f = Fs/N*(0:N/2);
figure;
plot(f,Pyy1(1:N/2+1))
hold on
plot(f,Pyy(1:N/2+1),'+r')
hold off
title(['Power spectral density. main freq - ',num2str(f_data)])
xlabel('Frequency (Hz)')
figure;
plot([normax(Rx) normax(Rx1) normax(Ry)])
title('acorr comparsion');
legend('unfilt est', 'filt est', 'proc');
xlabel('Frequency (Hz)')