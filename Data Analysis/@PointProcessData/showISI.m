function showISI( obj, varargin )
%Shows the spike interval histogram and its exponential fit

num2disp = 6;
[fact] = factor(num2disp);
x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
numFigures = ceil(obj.numChannels/num2disp);
cellIndx = (1:obj.numChannels);
for k = 1:numFigures
    figure(k);
    for n = 1:num2disp    
        if n+(k-1)*num2disp > obj.numChannels
            break;
        else
            dataNum = cellIndx(n+(k-1)*num2disp);
            try
                s(n) = subplot(x,y,n); 
                        intervals = diff(obj.spiketimes{dataNum});
                        nbins = ceil( 1.87 * (length(intervals))^0.4 ); % according to Michael Conn Methods in Neuroscience p.233
                        histfit(diff(obj.spiketimes{k}),nbins,'exponential');
                        mu = mean(intervals);
                        sigma = std(intervals);
                        title(s(n),['cell ',num2str(dataNum),' CV ',num2str(sigma/mu)]);            
            catch err
                break;
            end
        end 
        xlabel('t [sec]');ylabel('Spike counts');
    end
    legend('Data','Exponential - fit');
end

