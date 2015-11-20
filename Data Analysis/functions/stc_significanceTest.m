function [  ] = stc_significanceTest( rawStimuli , CP , var_signif, alpha, numIter )
%STC_significanceTest preforms significance test to STC as described in
%Schwartz, Pillow et al. 2006
% target - test if the importance of the significant kernels are greater
% than wht would expected by chance.
%         input: alpha - confidence interval (CI)
%                numIter - number of simulation for single CI estimation
%                rawStimuli - the stimulus matrix as described and calculated in compute_stc

% create Shift matrix for the CP shift - TODO CORRECT SHIFT MATRIX
mid = ceil( size(rawStimuli,1)/2 );
shiftMatrix = fliplr( [zeros( size(rawStimuli,1) - mid,size(rawStimuli,2)) ; eye( mid,size(rawStimuli,2) )] );
shiftfun = @(A,B) conv(A,B,'same');
CPshifted = cell2mat( bsxfun(@(x, y) shiftfun(x,y),CP, shiftMatrix) );

% iterative Monte Carlo simulation of CP shift and STE extraction
% for k = 1:numIter
%     CP_shift = 
%     X_squeeze = rawStimuli( logical(CP_shift),:);
%     CP_squeeze = CP( logical(CP) );
%     X = bsxfun(@times, X_squeeze ,CP_squeeze);
% end

end

