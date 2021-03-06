function [ ] = compareInNlPlane( obj )
% this function visually compares the derived lambdaCorr_forw from
% PreProcess to the analytical model derived in "the whitening effect in
% the nl plane.doc".

% extract notations from data structures
s = obj.stimulus.Yuncorr; 
h = {obj.STA.STA};

% prior of optimal STA gain
try
    N = load('staGain.mat');
    N = N.N;
catch exception
    N = ones(obj.numChannels,1);
end
h = arrayfun(@(k) h{k} * N(k),(1:obj.numChannels),'UniformOutput',false);

x = obj.CGP;
% theta = {obj.NlinKernel.estimation.fitParam};
nl = {obj.NlinKernel.estimation};
theta_ALL = cellfun(@(s) s.fitParam, nl,'UniformOutput', false);
R_lambd_b = obj.acorr.lambda_Corr_RAW;
[ R_lambd_f,R_lambd_b_BIN ] = deal(cell(numel(R_lambd_b),1));
maxlags = obj.maxlags;
dt = obj.dt;

% for each process compare
for k = 1:obj.numChannels
    
    % find statistical properties and compare assuming WSS
%     h{k} = h{k}; % the linear kernel should have a negative sign (?ASK SHY?)
    xk = conv(s, flipud(h{k}),'full' ); % x{k} is biased due to the 'same' property (the signal is not really WSS)
    xk = xk(1:length(s));
    E_xk = mean( xk ); 
    E_xk_est = mean( s ) * sum( h{k} );
    var_xk = var( xk );
    resampFactor = obj.stimulus.dt/(obj.dt);
    [R_s, lags] = xcorr(s,obj.maxlags/resampFactor,'unbiased');
    if size(h{k},1) == 1 % if h is not column vect -> transpose
        h{k} = h{k}';
    end
    kernHalfSize = round(length( h{k} )/2);
    conv_full = conv( flipud(h{k}),R_s,'full' );
    E_x_ttau = conv( conv_full(kernHalfSize:(end-kernHalfSize)) , h{k},'same' );
    var_xk_est = E_x_ttau(obj.maxlags/resampFactor+1) - E_xk_est.^2; 
    E_x_ttau = xcorr(conv( R_s,flipud(h{k}),'same'),obj.maxlags/resampFactor,'unbiased' );
    
    % downsampling to mutual low sampling rate
    zeroLag_index = obj.maxlags+1;
    R_lambd_b{k}(zeroLag_index) = [];
    lhs_R_lambd_b = binn( R_lambd_b{k}(zeroLag_index:end), resampFactor, @mean );
    R_lambd_b{k} = [flipud(lhs_R_lambd_b); lhs_R_lambd_b];
    maxlags = length(lhs_R_lambd_b);
    dt = obj.stimulus.dt;
    timeaxis = dt * (-maxlags:maxlags);
    timeaxis(maxlags+1) = []; 
    
    % find R_lambd_f according to formula
%     theta = fliplr(obj.NlinKernel.estimation(k).polyfit);
    theta = fliplr(theta_ALL{k});
    R_lambd_f{k} = theta(1)*theta(1) + 2*theta(1)*theta(2)*E_xk + 2*theta(1)*theta(3)*var_xk...
                + theta(2)*theta(2) .* E_x_ttau+ theta(3)*theta(3)... 
                .* (var_xk.^2 + 2.*E_x_ttau.^2) + 2*theta(2)*theta(3).*var_xk.*E_xk;  
    R_lambd_f{k}(maxlags + 1) = [];
    timeaxis = dt * (-maxlags:maxlags);
    timeaxis(maxlags+1) = [];    
    
%     % interpolation to mutual sampling rate
%     sampPoints =  lags .* obj.stimulus.dt;
%     queryPoints = linspace(sampPoints(1),sampPoints(end),...
%         (length(sampPoints)-1) * resampFactor + 1);
%     E_x_ttau_intrp = interp1( sampPoints, E_x_ttau ,queryPoints, 'spline' );
%     var_xk_est = E_x_ttau_intrp(obj.maxlags+1) - E_xk_est.^2; 
    
% %     find R_lambd_f according to formula
% %     theta = fliplr(obj.NlinKernel.estimation(k).polyfit);
%     theta = fliplr(theta_ALL{k});
%     R_lambd_f{k} = theta(1)*theta(1) + 2*theta(1)*theta(2)*E_xk + 2*theta(1)*theta(3)*var_xk...
%                 + theta(2)*theta(2) .* E_x_ttau_intrp + theta(3)*theta(3)... 
%                 .* (var_xk.^2 + 2.*E_x_ttau_intrp.^2) + 2*theta(2)*theta(3).*var_xk.*E_xk;  
    
    % remove zero lag   value for visual comparison
%     zeroLag_index = obj.maxlags+1;
%     R_lambd_b{k}(zeroLag_index) = [];
%     R_lambd_f{k}(zeroLag_index) = [];
%     timeaxis = dt * (-maxlags:maxlags);
%     timeaxis(maxlags+1) = [];    

%     find R_lambd_f numerically
    theta = theta_ALL{k};
    switch obj.NlinKernel(1).type      
        case 'poly' % poly case
            R_lambd_f{k} = xcorr( polyval(theta,xk),...
                maxlags,'unbiased'); % this is the crude numerical acorr with original h and s
        case 'exp'% exp case
            R_lambd_f{k} = xcorr( theta(1) * exp( theta(2) .* xk )...
            , maxlags,'unbiased'); % this is the crude numerical acorr with original h and s
    end
    
    R_lambd_f{k}(maxlags + 1) = [];
    
%     %** optional Binning
%         % binn for more poissonian like behaviour
%     binSize = 0.005; % possible values - 0.01 ~ 100 Hz 0.03
%     numBins = ceil(obj.maxlags*(obj.dt/binSize)); % define number of bins
%     binEdges = linspace( timeaxis(maxlags + 1), timeaxis(end),numBins);
%     [~,whichBin] = histc(timeaxis(maxlags + 1:end), binEdges);
%     R_right_side = R_lambd_b{k}(maxlags + 1:end);
%     for n = 1:numBins
%         flagBinMembers = (whichBin == n);
%         binMembers     = R_right_side(flagBinMembers);
%         R_BIN(n)       = mean(binMembers);
%     end
%     R_lambd_b_BIN{k} = interp1( [-fliplr(binEdges) binEdges], [fliplr(R_BIN) R_BIN] ,timeaxis, 'spline' );
end

% compare visually
plotProperties = struct('cellarray2',R_lambd_f','time',timeaxis,'method',...
    6,'num2disp',6,'title','R_{dN}(\tau) cell ','xlabel','\tau[sec]','ylabel'...
        ,'R(\tau) standertized','axis',[-0.5 0.5 -Inf Inf],'legendA',...
        'Data','legendB','LNP(STA) ');
plotCompare( R_lambd_b, plotProperties );
% plotCompare( R_lambd_b_BIN, plotProperties );

% plotCompare( R_lambd_b, R_lambd_f, timeaxis, 6, 'cell', normMethod, (1:cellstoshow*plotsPerFig));

% for figure in onenote
% plotCompare( R_lambd_b, R_lambd_f, timeaxis, 4, 'cell', 6, [4 5 7 8 15 20 21 23]); % figure 1
% plotCompare( R_lambd_b_BIN, R_lambd_f, timeaxis, 4, 'cell', 6, [4 5 7 8 15 20 21 23]); % figure 4
% plotCompare( R_lambd_b, R_lambd_f, timeaxis, 4, 'cell', 6, [4 5 7 8 15 20 21 23]); % figure 5

end
