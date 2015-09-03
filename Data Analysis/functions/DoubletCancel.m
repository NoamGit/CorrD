function [obj] = DoubletCancel( obj, varargin )
%Filters spikes according to empirical findings of refractory period and
%recovery function. for default function initialization see Berry & Meister
%1998

corrNoZero = ppcorr( obj, 'unbiased', 2, 'auto' );
for k = 1:obj.numChannels 
    % finding suspicious lag
    [~,IDX] = max( corrNoZero{k} );
    lags = (-obj.maxlags:obj.maxlags) .* obj.dt;
    lags(obj.maxlags) = [];
    doubletVal = abs(lags(IDX));
    doubletIndex = abs( diff( obj.spiketimes{k} ) -  obj.dt * floor(doubletVal/obj.dt)  ) < obj.dt;
    
    % filter conseq doublets
    y = zeros(size(doubletIndex));
    for n = 2:length(doubletIndex)
        y(n) = logical(doubletIndex(n) * ( (y(n-1)+doubletIndex(n))/2 - 1 ));
    end
    doubletIndex = [ false; logical(y) ];
    obj.spiketimes{k}(doubletIndex) = [];
end
end

