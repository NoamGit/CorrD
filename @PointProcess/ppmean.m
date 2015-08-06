function [mfr] =ppmean(obj)

% Ex=ppmean(obj.spiketimes) estimates the mean rates of the point processes
% Ex - a scalar or a cell array of estimated mesn rates
% obj.spiketimes - a vector or a cell array of the event obj.spiketimes of the processes

vFlag=false;
if ~iscell(obj.spiketimes)
    vFlag=true;
    % reshape obj.spiketimes vector to be a column vector 
    % (also returns an error if obj.spiketimes is not a 1-D vector)
    obj.spiketimes=reshape(obj.spiketimes, length(obj.spiketimes), 1);
    obj.spiketimes=mat2cell(obj.spiketimes, size(obj.spiketimes, 1), size(obj.spiketimes, 2));
end

mfr=cellfun(@(X) length(X)/(max(X)),obj.spiketimes);

% return a scalar (and not a cell array) if obj.spiketimes was a vector
if vFlag
    mfr=mfr{1};
end
