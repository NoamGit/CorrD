X = (0:1e-3:2);
m1 = 0.3; s1 = 0.2; n1 = 1;
Y1 = n1 .* normcdf(X, m1,s1);
m2 = -8; s2 = 4; n2 = 3;
Y2 = n2.*exp(m2+s2.*X);
m3 = 1.5 ;s3 = 0.1; n3 = 7;
Y3 = n3 .* normcdf(X, m3,s3);
plot(X,Y1+Y3,'--d',X,Y1,'k',X,Y2,'--r',X,Y3,'m');

%%
raw_corr_comparison;

pnum = 19; % process num
maxlags = process.maxlags;
timeaxis = fliplr(process.dt .* (1:maxlags));
Y = flipud(process.acorr.lambda_Corr_RAW{pnum}(maxlags+2:end)); % lambda corr
stim = resample( process.stimulus.Yuncorr, process.Fs, Fs_STIM );
sta = resample( process.STA(pnum).realSTA, process.Fs, Fs_STIM );
X = conv(sta, stim); % CGP
X_corr_forw = xcov(X, maxlags, 'unbiased');
X_corr_forw = flipud(X_corr(maxlags+2:end));

%%
mu = 50 ; sigma =3;
X_corr_exp = (1/(sigma^2)) * log( Y/0.01); % lambda corr est

n1 = 0.08 * sigma; mu1 = 70; sig1 = 3 * sigma;
pd = makedist('Normal',mu1,sig1);
X_corr_back = (1/(sigma^2)) * log( Y/0.01 ) + n1 * normcdf( Y, mu1, sig1 ) ; 
% X_corr_back = (n3 .* icdf(pd, 0.3* Y/(max(Y))) ) ; 

figure(1);
plot(timeaxis, normax(X_corr_back-mean(X_corr_back)), timeaxis,...
    normax(X_corr_forw-mean(X_corr_forw)),timeaxis,normax(X_corr_exp-mean(X_corr_exp)));
legend('PROBIT CGP back acorr','generated CGP forward acorr',' EXP CGP Back acorr')
axis([timeaxis(end) timeaxis(1) -2 2]);

figure(2);
plot(timeaxis,(X_corr_back-mean(X_corr_back)));
legend('PROBIT CGP back acorr')

figure(4)
plot((1/(sigma^2)) * log( Y/1e2 ) - n1 * normpdf( Y, mu1, sig1 ))
