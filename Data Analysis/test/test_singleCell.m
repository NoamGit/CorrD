%% case study single cell

k = 11;
s = process.stimulus.Yraw;
s = s-mean(s);
h = process.STA(k).STA;
x = process.CGP{k};
theta = process.NlinKernel(k).estimation.fitParam;
lags = process.acorr.lags;
R_dN = process.acorr.lambda_Corr_RAW{k};
if(strcmp( process.NlinKernel(1).type,'exp'))
    nl = @(z) theta(1) * exp( theta(2) .* z );
elseif(strcmp( process.NlinKernel(1).type,'poly'))
    nl = @(z) polyval(theta,z);
else
    disp('nl type isn''t defined');
end
[R_s, lags] = xcorr(s,process.maxlags/resampFactor,'unbiased');
E_x_ttau = conv( conv( R_s,flipud(h),'same' ), h,'same' );
resampFactor = process.stimulus.dt/(process.dt);
zeroLag_index = process.maxlags+1;
R_dN(zeroLag_index) = [];
lhs_R_dN = binn( R_dN(zeroLag_index:end), resampFactor, @mean );
R_dN = [flipud(lhs_R_dN); lhs_R_dN];
R_LNP = xcorr( nl(x), process.maxlags,'unbiased');
R_LNP(zeroLag_index) = [];

maxlags = length(lhs_R_LNP);
dt = process.stimulus.dt;
timeaxis = dt * (-maxlags:maxlags);
timeaxis(maxlags+1) = []; 
E_x_ttau(maxlags+1) = [];

figure; 
subplot(221) 
plot(timeaxis, normalize(R_LNP_resap,6));
axis([-0.5 0.5 -Inf Inf]);
title('R _{LNP}')
subplot(222)
plot(timeaxis,normalize(E_x_ttau,6));
title('E[x_{t}x_{t+1}]')
axis([-0.5 0.5 -Inf Inf]);
subplot(223)
plot(timeaxis,normalize(R_dN,6));
title('R _{dN}')
axis([-0.5 0.5 -Inf Inf]);