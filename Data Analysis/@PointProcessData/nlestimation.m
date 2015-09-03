function obj = nlestimation( obj,varargin ) 
% nlestimation() estimates the non linearity out of the spike train. CGP is
% calculated according to the stimulus and the STA.
% the STA must be obtained before running this function (calcSTA());
% stimulus
%       inputs:  varargin{1}  - polyorder  for fitting the obtained data   

showFlag = 0;
polyorder = 0;
resampleFlag = 0;
if nargin > 0
    polyorder = varargin{1};
    showFlag = varargin{2};
    resampleFlag = any(strcmp(varargin,'resample stimulus')); 
end

CP = obj.CP;
nl_est = cell(obj.numChannels,2);
[ image_bin, nl_bin, param ] = deal( cell(obj.numChannels,1) );
for k = 1:obj.numChannels
        if ~resampleFlag % if stimulus was not resampled interpolate CGP
            stim = [obj.stimulus.Yuncorr;  obj.stimulus.Yraw(length(obj.stimulus.Yuncorr)+1)-mean(obj.stimulus.Yraw)];
            sampVal = conv( stim, obj.STA(k).realSTA, 'same');
            sampPoints =  (0:obj.stimulus.dt:obj.T);
            queryPoints = (0:obj.dt:obj.T);
            CGP = interp1( sampPoints, sampVal ,queryPoints, 'spline' );
            obj.CGP{k} = CGP(1:end-1)'; % removing interp bias
            CP = obj.CP;
        else % everything is already resampled
            obj.CGP{k} = conv( obj.stimulus.Yuncorr, obj.STA(k).realSTA, 'same');
%             time =  (0:obj.stimulus.dt:obj.T-obj.stimulus.dt); % this can
%             be used for no interp nor resampling nl estimation
%             CP{k} =  histc(obj.spiketimes{k}, time);
        end
        
        % normalize the linear prediction and weighted binning of the functions domain
        u = obj.CGP{k}/std(obj.CGP{k}); 
        p = prctile(u',linspace(0,100,101)); 
        p(end)=p(end)+eps;

        % plug CGP values to binned domain and keep the order
        [N,bin]=histc(u,p); 
        N=N(1:end-1);

        % unwanted bin index
        unwant = find(bin > length(N)); 
        u(unwant) = [];
        CP{k}(unwant) = [];
        bin(unwant) = [];
        
        domain_bin = accumarray(bin,u,[],@mean); % the domain (X axis) of the non linearity
        image_bin{k} = accumarray(bin,CP{k},[],@mean); % the image (Y axis) find how many spikes for each bin
        param{k} = polyfit(domain_bin,image_bin{k},polyorder); % calculate polynomial fit
        nl_est{k,1} = domain_bin;
        nl_est{k,2} = polyval(param{k},domain_bin);
        nl_bin{k} = bin;
end

obj.NlinKernel.estimation = struct('domain',nl_est(:,1),'image',nl_est(:,2),'bin', nl_bin,'polyfit',param(:));

if showFlag
num2disp = 4;
[fact] = factor(num2disp);
x = prod(fact(1:ceil(length(fact)/2))); y =  prod(fact(ceil(length(fact)/2)+1:end));
numFigures = ceil(obj.numChannels/num2disp);
    for m = 1:numFigures
        figure();
        for n = 1:num2disp
            dataNum = (n+(m-1)*num2disp);
            s(n) = subplot(x,y,n);
            try
                plot(nl_est{dataNum,1},nl_est{dataNum,2},nl_est{dataNum,1},image_bin{dataNum,:}, '+r');
                title(s(n),['nl number ', num2str(dataNum)]);
            catch err
            end
        end
    end
end
end

    
