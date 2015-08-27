function [ ] = compare_whitening_effect_1( obj )
% this function visually compares the derived lambdaCorr_forw from
% PreProcess to the analytical model derived in "the whitening effect in
% the nl plane.doc"

% extract notations from data structures
s = obj.stimulus.Yuncorr; 
h = {obj.STA.realSTA};
x = obj.CGP;
theta = {obj.NlinKernel.estimation.polyfit};
R_lambd_b = obj.acorr.lambda_Corr_RAW;
R_lambd_f = cell(numel(R_lambd_b),1);

% for each process compare
for k = 1:obj.numChannels
    
    % find statistical properties and compare assuming WSS
    xk = conv( h{k}, s, 'full' ); % x{k} is biased due to the 'same' property (the signal is not really WSS)
    E_xk = mean( xk ); 
    E_xk_est = mean( s ) * sum( h{k} );
    var_xk = var( xk );
    resampFactor = obj.stimulus.dt/(obj.dt);
    [R_s, lags] = xcorr(s,obj.maxlags/resampFactor,'unbiased');
    E_x_ttau = conv( conv( R_s,flipud(h{k}),'same' ), h{k},'same' );
    
    % interpolation to mutual sampling rate
    sampPoints =  lags .* obj.stimulus.dt;
    queryPoints = linspace(sampPoints(1),sampPoints(end),...
        length(sampPoints) * resampFactor);
    E_x_ttau_intrp = interp1( sampPoints, E_x_ttau ,queryPoints, 'spline' );
    var_xk_est = E_x_ttau_intrp(obj.maxlags+1) - E_xk_est.^2; 
    
%     find R_lambd_f according to formula
    theta = fliplr(obj.NlinKernel.estimation(k).polyfit);
    R_lambd_f{k} = theta(1)*theta(1) + 2*theta(1)*theta(2)*E_xk + 2*theta(1)*theta(3)*var_xk...
                + theta(2)*theta(2) .* E_x_ttau_intrp + theta(3)*theta(3)... 
                .* (var_xk.^2 + 2.*E_x_ttau_intrp.^2) + 2*theta(2)*theta(3).*var_xk.*E_xk;  
    
%     R_lambd_f{k} = xcorr( polyval(obj.NlinKernel.estimation(k).polyfit,...
%     xk), obj.maxlags,'unbiased'); % this is the crude numerical acorr
%     with original h and s
end

% compare visually
plotCompare( R_lambd_b, R_lambd_f, obj.acorr.lags, 6, 'cell', 1, (1:24));
 
end
