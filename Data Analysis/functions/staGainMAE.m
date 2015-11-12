function [ cost ] = staGainMAE( gain, s, h, nl, maxlags, R_dN )
%staGainObjfun - MAE measurment for R_dN and R_LNP differenc efor different
%STA gain. used in test_single cell for optimal gain retrival

x = conv(s, flipud(gain .* h),'full' ); 
x = x(1:length(s));

R_LNP = xcorr( nl(x), maxlags,'unbiased');
R_LNP(maxlags+1) = [];

cost = sum( abs( normalize(R_dN,6) - normalize(R_LNP,6) ) );
end

