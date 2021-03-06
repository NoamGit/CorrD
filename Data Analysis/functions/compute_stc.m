function [ eig_vec, mu, explained, eig_val, X, rawStimuli ] = compute_stc( stimulus, CP, numSamples, sta, method)
%compute_stc Calculates the spike-triggered covariance for a neuron that
%           is driven by a stimulus defined in stim. The spike-
%           inputs: 
%                   rwe - window's right edge 
%                   method -  of sta axis is suppressed with 'suppress sta'
%                   we project out the sta axis from the stimulus ensemble
%                   subspace
%  Description:     
%           1. finds X - the 'Spike Triggered Ensemble' ('STE' - a collection of all
%           stimulus paths that happen before a spike). 
%           2. calculates X * X' the covariance matrix of the STE
%  ref - Spiked neural characterization - Shwartz, pillow et al 2006
    
    % 1. find STE
    rawStimuli = repmat(stimulus,1,numSamples);
    mid = ceil( size(rawStimuli,1)/2 );
    shiftMatrix = fliplr( [zeros( size(rawStimuli,1) - mid,numSamples) ; eye( mid,numSamples  )] );
    shiftfun = @(A,B) conv(A,B,'same');
    rawStimuli = cell2mat( arrayfun(@(k) shiftfun(rawStimuli(:,k),shiftMatrix(:,k)),...
        (1:numSamples),'UniformOutput',false));
    X_squeeze = rawStimuli( logical(CP),:);
    CP_squeeze = CP( logical(CP) );
    X = bsxfun(@times, X_squeeze ,CP_squeeze);
    if strcmp(method,'suppress sta')
        if (size(sta,1) == 1)
            sta = sta';
        end
        X_StaOutProjected = X - bsxfun(@times,sta', X * sta)./(norm(sta).^2);
    else
        X_StaOutProjected = X;
    end
    
    % 2. cacl COV matrix
    [eig_vec, score, eig_val, tsqr ,explained, mu ] = pca(X_StaOutProjected);
    
    % omitting out the last value if supress sta is in use
    if strcmp(method,'suppress sta')
        eig_vec = eig_vec(:, 1:end-1);
        eig_val = eig_val(1:end-1);
        explained = explained(1:end-1);
    end
    
%     disp(['explained 1-10 ', mat2str(explained(1:10)')]);
end



