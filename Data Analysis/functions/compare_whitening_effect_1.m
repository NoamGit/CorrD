function [ ] = compare_whitening_effect_1( obj )
% this function visually compares the derived lambdaCorr_forw from
% PreProcess to the analytical model derived in "the whitening effect in
% the nl plane.doc"

% extract notations from data structures
s = obj.stimulus.Yuncorr; 
h = {obj.STA.realSTA};
x = obj.CGP;
theta = {obj.NlinKernel.estimation.polyfit};
R_lambd_f = obj.acorr.R_lambd_f;

% for each process compare
for k = 1:obj.numChannels
    
    % find statistical properties and compare assuming WSS
    xk = conv( h{k}, s, 'full' ); % x{k} is biased due to the 'same' property (the signal is not really WSS)
    E_xk = mean( xk ); 
    E_xk_est = mean( s ) * sum( h{k} );
    var_xk = var( xk );
    var_xk_est = conv( h; 
    
end

end
