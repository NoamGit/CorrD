function [ ] = JointFilterResponseProb( k0, k1, nbins, X ,figureHandle )
%JointFilterResponseProb finds and plots the joint filter response
%probabilty matrix for 2 filters.
%         input: k0 - filter 1
%                k1 - filter 2
%                binSize - the quantization bin size of the responses
%                X - the stimulus ensemble matrix
% test with code :  binmap1 = [1 1 2 3 2 1 ]';dotK0_norm = [1 1 -1 0 0 1
% ]'; dotK1_norm = [-1 -1 0 1 0 -1]';nbins = 3;

% find reponses to filters by simple dot product (like convolution)
if (size(k0,1) == 1) % flipping dimenstion
    dotK0 = X * k0';
    dotK1 = X * k1';
else
    dotK0 = X * k0;
    dotK1 = X * k1;
end

% normalize data to [-1,1]
dotK0_norm = 2 * norm_nc(dotK0,5) - 1;
dotK1_norm = 2 * norm_nc(dotK1,5) - 1;

% binn into nbins bins
[~, edges, ~] = histcounts(dotK0_norm, nbins);
[~, ~, binmap1] = histcounts(dotK1_norm, nbins);

% calculates the histogram conditions on each binmap0
[~,sortmap] = sort(dotK1_norm, 'ascend');
condition = binmap1(sortmap);
k0_sorted = dotK0_norm(sortmap);
jointMatrix = arrayfun(@(x) fliplr(hist( k0_sorted( condition == x ) ,edges ) ), (1:numel(edges)), 'UniformOutput',false);
jointMatrix = cell2mat(jointMatrix')';

% plot joint conditial  matrix
% Create axes
axes1 = axes('Parent',figureHandle,'Xlim',[1 40],'Ylim',[1 40],'YTickLabel',{'1','0','-1'},'YTick',[0 20 40],...
'XTickLabel',{'-1','0','1'},'XTick',[0 20 40],'Layer','top','YDir','reverse',...
'DataAspectRatio',[1 1 1]);
box(axes1,'off');
hold(axes1,'on');

% Create image
image(get(imagesc(jointMatrix),'CData'),'Parent',axes1,'CDataMapping','scaled');
end

