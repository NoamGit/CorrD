%% Report for AcorrSTA_Sim_2_1
% run first [est, proc] = AcorrSTA_Sim_2_1 for  n = 100 instances

[Rsquare, RMSE_proc] = deal( zeros(length(est), 1) );
[param_data, param_est] = deal( zeros(2,length(est)) );
for n = 1:length(est)
    RMSE_proc(n) = est{n}.PEresults.estimationRMS;
    Rsquare(n) = est{n}.PEresults.GoodnessStatistics.rsquare;
    param_est(:,n) = est{n}.PEresults.paramEst;
    param_data(:,n) = est{n}.PEresults.paramData; % [sig1, mu1, pi, f ]
end

n = length(est);
flag = 0;
badI = find(RMSE_proc > 0.5); % good estimation/data indexes
goodI = setdiff((1:n),badI); % bad estimation/data indexes
good_samp_idx = goodI(randperm(length(goodI),4));
good_samp = est(good_samp_idx);
if ~isempty(badI)
%     bad_samp_idx = badI(randperm(length(badI),2));
    bad_samp_idx = badI(1:2);
    flag = 1;
end
%% Display numerous examples

% Display good samples of the estimated acorr_CGP and the real data
cellfun(@(x,y) plotEstimation(x,y), est(good_samp_idx), proc(good_samp_idx) );

% Display bad samples of the estimated acorr_CGP and the real data
if flag
cellfun(@(x,y) plotEstimation(x,y), est(bad_samp_idx), proc(bad_samp_idx) );end;
%% Analysis of entire Data measurments

RMSE_bar = mean(RMSE_proc(goodI));
RMSE_std = std(RMSE_proc(goodI));
Rsquare_bar = mean(Rsquare(goodI));
Rsquare_std = std(Rsquare(goodI));
%% Histogram of R^2

hHist1 = figure;
set(hHist1, 'Position', [100, 150, 1100, 650]);
xvalues = 0:1;
hist(Rsquare(Rsquare >=  0),50,xvalues);
set(gca,'FontSize',14);
title('Histogram of \itR^2 values');
annotation(hHist1,'textbox',[0.2 0.78 0.2 0.1],... % textbox [ x y w h]
'String',{['\itR^2 \pm \it\sigma = ',num2str(Rsquare_bar),' + ',num2str(Rsquare_std)]},'FitBoxToText','on',...
'EdgeColor','none','FontSize',16);
hHist1 = findobj(gca,'Type','patch');
set(hHist1,'FaceColor',[0 .5 .5],'edgeColor','k');
%% Histogram of RMSE - 

% Note: The parameter difference is divided by the mean of the jointed
% parameter values. this division emphasise the similarity between the
% values and supresses the case where we have simply a small difference
% without similarity.
hHist2 = figure;
set(hHist2, 'Position', [100, 150, 1100, 650]);
hist(RMSE_proc,50);
set(gca,'FontSize',14);
title('Histogram of \itRMSE values');
annotation(hHist2,'textbox',[0.2 0.78 0.2 0.1],... % textbox
'String',{['\itRMSE \pm \it\sigma = ',num2str(RMSE_bar),' + ',num2str(RMSE_std)]},'FitBoxToText','on',...
'EdgeColor','none','FontSize',16);
hHist2 = findobj(gca,'Type','patch');
set(hHist2,'FaceColor',[0 .5 .5],'edgeColor','k');

