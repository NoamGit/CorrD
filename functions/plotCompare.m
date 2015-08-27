function [ ] = plotCompare(cellarray1, varargin)
% vizualizes graphically 2 cell arrays
% input: num2disp number of plots to display in a figure

if nargin>2
    cellarray2 = varargin{1};
    timeaxis = varargin{2};
    num2disp = varargin{3};
    title_arg = varargin{4};
    normFlag = varargin{5};
    cellIndx = varargin{6};
else
    timeaxis = (1:length(cellarray1{1}));
    num2disp = 4;
    title_arg = 'No title';
    normFlag = 1;
    cellIndx = (1:numel(cellarray1));
end

[fact] = factor(num2disp);
x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
numFigures = ceil(numel( cellarray1 )/num2disp);   
for k = 1:numFigures
    figure(k);
    for n = 1:num2disp
        dataNum = cellIndx(n+(k-1)*num2disp);
        try
            s(n) = subplot(x,y,n);
            if exist('cellarray2','var')
                if normFlag
                plot(timeaxis, normax(cellarray1{dataNum}-mean(cellarray1{dataNum}))...
                    , timeaxis, normax(cellarray2{dataNum}-mean(cellarray2{dataNum})),'r');
                title(s(n),[title_arg,' ',num2str(dataNum)]);
%                     plot(timeaxis, normax(cellarray1{dataNum})...
%                     ,'--', timeaxis, normax(cellarray2{dataNum}),'r');
                else
                    plot(timeaxis, cellarray1{dataNum},'--', timeaxis, cellarray2{dataNum},'r');
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                end
            else
                 if normFlag
    %             plot(timeaxis, normax(cellarray1{dataNum}-mean(cellarray1{dataNum}))...
    %                 ,'--g', timeaxis, normax(cellarray2{dataNum}-mean(cellarray2{dataNum})),'r');
                    plot(timeaxis, normax(cellarray1{dataNum}));
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                else
                    plot(timeaxis, cellarray1{dataNum},'--');
                    title(s(n),[title_arg,' ',num2str(dataNum)]);
                 end
            end
            
        catch err
            break;
        end
        axis([-0.3 0.3 -0.2 1]);
%         axis([timeaxis(1) timeaxis(end) -inf inf]);
%         title(s(n),[title,' ',num2str(dataNum)]); 
        xlabel('\tau [sec]');ylabel('corr');
    end
%     legend('CGP corr','STA corr');
end
end

