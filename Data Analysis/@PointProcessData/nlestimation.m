function [ obj, domain_bin] = nlestimation( obj,varargin ) 
% nlestimation() estimates the non linearity out of the spike train. CGP is
% calculated according to the stimulus and the STA.
% the STA must be obtained before running this function (calcSTA());
% stimulus
%       inputs:  varargin{1}  - polyorder  for fitting the obtained data   

showFlag = 0;
polyorder = 0;
resSimFlag = 0;
amp = 1;
syntDataFlag = 0;
if nargin == 3
    polyorder = varargin{1};
    showFlag = varargin{2};
    method = varargin{3};
    resSimFlag = any(strcmp(varargin,'resampled stimulus')); 
elseif nargin > 3
    polyorder = varargin{1};
    showFlag = varargin{2};
    amp = varargin{3};
    method = varargin{4};
    syntDataFlag = any(strcmp(varargin,'Synt Data'));
    resSimFlag = any(strcmp(varargin,'resampled stimulus'));
    resCPFlag = any(strcmp(varargin,'resample CP'));
end

CP = obj.CP;
nl_est = cell(obj.numChannels,2);
[ image_bin, nl_bin, param ] = deal( cell(obj.numChannels,1) );
for k = 1:obj.numChannels
        if ~resSimFlag && ~syntDataFlag && ~resCPFlag % case 1 - stimulus was not resampled and data isn't syntetic  -> interpolate CGP 
            stim = [obj.stimulus.Yuncorr;  obj.stimulus.Yraw(length(obj.stimulus.Yuncorr)+1)-mean(obj.stimulus.Yraw)];
            sampVal = conv( stim, flipud(obj.STA(k).STA), 'full');
            sampVal = sampVal(1:end-length(obj.STA(k).STA)+1);
            sampPoints =  (0:obj.stimulus.dt:obj.T);
            queryPoints = (0:obj.dt:obj.T);
            CGP = interp1( sampPoints, sampVal ,queryPoints, 'spline' );
            obj.CGP{k} = CGP(1:end-1)'; % removing interp bias
            CP = obj.CP;
            u = obj.CGP{k};
            %         u = obj.CGP{k}/std(obj.CGP{k}); 
        elseif resCPFlag && ~resSimFlag && ~syntDataFlag
            stim = [obj.stimulus.Yuncorr;  obj.stimulus.Yraw(length(obj.stimulus.Yuncorr)+1)-mean(obj.stimulus.Yraw)];
            CGP = conv( stim, flipud(obj.STA(k).STA), 'full');
            CGP = CGP(1:end-length(obj.STA(k).STA)+1);
            binningFactor = length(obj.CP{k})/(numel(CGP)-1);
            cp_binn = binn(obj.CP{k}, binningFactor, @sum);
            CP{k} = [cp_binn ; 0];
            u = CGP;
%             u = CGP/std(CGP);
        elseif resSimFlag && ~syntDataFlag % case 2 - everything is already resampled and not syntetic data
            cgp =  conv( obj.stimulus.Yuncorr, amp * flipud(obj.STA(k).STA), 'full');
            obj.CGP{k} = cgp(1:end-length(obj.STA(k).STA)+1);
            u = obj.CGP{k};
            %         u = obj.CGP{k}/std(obj.CGP{k}); 
%             time =  (0:obj.stimulus.dt:obj.T-obj.stimulus.dt); % this can
%             be used for no interp nor resampling nl estimation
%             CP{k} =  histc(obj.spiketimes{k}, time);   
        else 
            cgp =  conv( obj.stimulus.Yuncorr, amp * flipud(obj.STA(k).STA), 'full');
            obj.CGP{k} = cgp(1:end-length(obj.STA(k).STA)+1);
            u = obj.CGP{k};
        end
        
        % normalize the linear prediction and weighted binning of the functions domain
%         u = obj.CGP{k}/std(obj.CGP{k}); 
        p = prctile(u',linspace(0,100,101)); 
        p(end)=p(end)+eps;

        % plug CGP values to binned domain and keep the order. bin is a
        % mapping of the u values to ther percentile's centers ( it's hist
        % should be approximatly uniform)
        [N,bin]=histc(u,p); 
        N=N(1:end-1);

        % unwanted bin index
        unwant = find(bin > length(N)); 
        u(unwant) = [];
        CP{k}(unwant) = [];
        bin(unwant) = [];
        
        domain_bin = accumarray(bin,u,[],@mean); % the domain (X axis) of the non linearity
        image_bin{k} = accumarray(bin,CP{k},[],@mean); % the image (Y axis) find how many spikes for each bin
%         scatter(domain_bin,image_bin{k});
        
        switch method
            case 'poly'
                param{k} = polyfit(domain_bin,image_bin{k},polyorder); % calculate polynomial fit
                nl_est{k,1} = domain_bin;
                nl_est{k,2} = polyval(param{k},domain_bin);
                nl_bin{k} = bin;
                
            case 'exp'
                expfit_OBJ = fit(domain_bin,image_bin{k},'exp1');
                param{k} = [expfit_OBJ.a expfit_OBJ.b];
                nl_est{k,1} = domain_bin;
                nl_est{k,2} = expfit_OBJ.a * exp( expfit_OBJ.b .* domain_bin ); 
                nl_bin{k} = bin;
        end
        
        % CAUTION - might not work with real (not synt) data
        obj.NlinKernel(k).estimation = struct('domain',nl_est(k,1),'image',nl_est(k,2),'bin', nl_bin(k),'fitParam',param(k));
end

if showFlag
num2disp = 6;
[fact] = factor(num2disp);
x = prod(fact(1:ceil(length(fact)/2))); 
y =  prod(fact(ceil(length(fact)/2)+1:end));
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

    
