%% case study single cell
% exclude cells 1 7 9 10 13 14 16 21 22 

N = ones(23,1);
for k = (1:23)
    s = process.stimulus.Yuncorr;
    s = s-mean(s);
    h = 1*process.STA(k).STA;
    % x = process.CGP{k};
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
    resampFactor = process.stimulus.dt/(process.dt);

    x = conv(s, flipud(h),'full' ); % x{k} is biased due to the 'same' property (the signal is not really WSS)
    x = x(1:length(s));
    [R_s, lags] = xcorr(s,process.maxlags/resampFactor,'unbiased');
    E_x_ttau = conv( conv( R_s,flipud(h),'same' ), h,'same' );
    zeroLag_index = process.maxlags+1;
    R_dN(zeroLag_index) = [];
    lhs_R_dN = binn( R_dN(zeroLag_index:end), resampFactor, @mean );
    R_dN = [flipud(lhs_R_dN); lhs_R_dN];
    maxlags = length(lhs_R_dN);

    R_LNP = xcorr( nl(x), maxlags,'unbiased');
    R_LNP(maxlags+1) = [];

    dt = process.stimulus.dt;
    timeaxis = dt * (-maxlags:maxlags);
    timeaxis(maxlags+1) = []; 
    E_x_ttau(maxlags+1) = [];
    %% optimize STA gain

    x0 = 1;
    objfun = @(N) staGainMAE(N, s, h, nl, maxlags, R_dN );
    N( k ) = fminsearch( objfun ,x0 );

    x = conv(s, flipud(N( k ) .* h),'full' ); 
    x = x(1:length(s));
    R_LNP = xcorr( nl(x), maxlags,'unbiased');
    R_LNP(maxlags+1) = [];
    %% plot

    figure(1); 
    subplot(221) 
    plot(timeaxis, normalize(R_LNP,6));
    axis([-0.8 0.8 -Inf Inf]);
    title('R _{LNP}')
    subplot(222)
    plot(timeaxis,normalize(E_x_ttau,6));
    title('E[x_{t}x_{t+1}]')
    axis([-0.8 0.8 -Inf Inf]);
    subplot(223)
    plot(timeaxis,normalize(R_dN,6));
    title('R _{dN}')
    axis([-0.8 0.8 -Inf Inf]);
    subplot(224)
    plot(timeaxis,normalize(R_dN,6),timeaxis, normalize(R_LNP,6),'--r');
    title('R _{dN}');
    axis([-0.8 0.8 -Inf Inf]);
end
