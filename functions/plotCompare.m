function [ ] = plotCompare(cellarray1, varargin)
% vizualizes graphically 2 cell arrays
% input: num2disp number of plots to display in a figure
% for normalization method (varargin{5}) see normalize.m

timeaxis = (1:length(cellarray1{1}));
num2disp = 4;
title_arg = 'No title';
method = 0;
cellIndx = (1:num2disp * ceil(numel(cellarray1)/num2disp));

if nargin == 2
    timeaxis = varargin{1};
elseif nargin == 3
    timeaxis = varargin{1};
    cellIndx = varargin{2};
elseif nargin > 2
    cellarray2 = varargin{1};
    timeaxis = varargin{2};
    num2disp = varargin{3};
    title_arg = varargin{4};
    method = varargin{5};
    cellIndx = varargin{6};
end

if ~iscell(timeaxis)
    timeRepMat = repmat(timeaxis,numel(cellarray1),1);
    timeaxis = mat2cell(timeRepMat, ones(numel(cellarray1),1),size(timeaxis,2));
end

[fact] = factor(num2disp);
x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
numFigures = ceil(numel( cellIndx )/num2disp);   
for k = 1:numFigures
%     figure();
    figure(k);
    for n = 1:num2disp
        try
            dataNum = cellIndx(n+(k-1)*num2disp);
            s(n) = subplot(y,x,n);
            if exist('cellarray2','var') % plot - 2 signals to compare
                if method
                plot(timeaxis{dataNum}, normalize(cellarray1{dataNum},method)...
                    , timeaxis{dataNum}, normalize(cellarray2{dataNum},method),'.r');
                title(s(n),[title_arg,' ',num2str(dataNum)]);
                else
                    plot(timeaxis{dataNum}, cellarray1{dataNum},'--', timeaxis{dataNum}, cellarray2{dataNum},'r');
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                end
            else % plot - only one signal 
                 if method
                    plot(timeaxis{dataNum}, normalize(cellarray1{dataNum},method));
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                else
                    plot(timeaxis{dataNum}, cellarray1{dataNum});
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                 end
            end
            
        catch err
            break;
        end
%         axis([-.5 .5 -Inf Inf]);
        axis([timeaxis{dataNum}(1) timeaxis{dataNum}(end) -inf inf]);
%         title(s(n),[title,' ',num2str(dataNum)]); 
        xlabel('Linear prediction');ylabel('Expected firing rate (Hz)');
    end
%     legend('CGP corr','STA corr');
%     legend('R_{dN}{(\tau)}','R_{\lambda}{(\tau)}')
end
end

