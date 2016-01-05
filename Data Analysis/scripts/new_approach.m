% Script - new approach 2/1/16

% load data
try
    data = load('C:\Users\noambox\Documents\CorrD\SourceData\Ronen\SpikeTimeGauss1B.mat'); % change path  
catch exception
    try
         data = load('C:\Users\Noam\Documents\GitHub\CorrD\SourceData\Ronen\SpikeTimeGauss1B.mat'); % change path
    catch exception
         disp(exception.message);
        display('file was not found. Please load manually...');
        [fileN, pathN ] = uigetfile('*.mat','Select Data file');    % C:\Users\Noam\Documents\GitHub\...SpikeTimeGauss1B.mat
        data = load([ pathN, fileN ]); % change path
    end
end 

% notice - max spike time 947 sec
% clean & arrange stimulus
s_dt = 30;
s_frame = data.W;
s = s_frame(1:947*s_dt,:);
flag = find(ismember(s(:,1),[175 255]));
flag_times = flag / s_dt; 
flag = flag( diff([flag(1); flag]) > 1);
experiment_samples = [[1; flag(1:end-1)] flag+1];
w = s(:,3);
clear s;
s = arrayfun(@(x,y) w(x:y), experiment_samples(:,1),experiment_samples(:,2),'Uniformoutput',false);
nTrails = numel(s);

% clean ang arrange counting process
% r_all is the ensemble of all cell reactions to a part of the stimulus
r_dt = 1e4;
experiment_samples_r = (experiment_samples/s_dt) * r_dt ;
r = {data.TT.sp};
r = cellfun(@(x) (arrayfun(@(y,z) x(y < x & x < z) ,...
    experiment_samples_r(:,1),experiment_samples_r(:,2),...
    'Uniformoutput',false)),r,'Uniformoutput',false);
nCell = numel(r);
r_all = [r{:}];
r_all = r_all(:);
r_all = arrayfun(@(x) r_all(x:numel(s):end), (1:numel(s)),'UniformOutput', false)';
[T,N,t] = deal(zeros(nTrails,1));
clear r;
for k = 1:nTrails
    t_start = (k-1) * length(s{k})*1/s_dt;
    t_end = k * length(s{k})*1/s_dt; % end of stimulus
    N = ceil( t_end * 450 );  
    t = linspace(t_start,t_end,N);
    r{k} = cellfun(@(x) histc(x./r_dt, t), r_all{k} ,'UniformOutput',false);
end
r = r';
save('.\SourceData\Ronen\ronens_data.mat','s','r')