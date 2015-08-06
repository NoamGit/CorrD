%% Report for estSTA_Sim_2_2
% run first [ est ] = estSTA_Sim_2_2 for  n = 100 instances

load('D:\# Projects (Noam)\# CorrD\code\trunk\Workspace\WS_estSTA_Sim_2_2.mat');

[Rsquare, RMSE_proc] = deal( zeros(length(est), 1) );
[param_data, param_est] = deal( zeros(4,length(est)) );
for n = 1:length(est)
    Rsquare(n) = est(n).PEresults.GoodnessStatistics.rsquare;
    param_est(:,n) = est(n).PEresults.paramEst;
    param_data(:,n) = est(n).PEresults.paramData; % [sig1, mu1, pi, f ]
    SSE_param(n) = sum( (param_est(:,n)- param_data(:,n)).^2 );
end

n = length(est);
flag = 0;
num2disp = 4;
badI = find(Rsquare < 0.95); % good estimation/data indexes
goodI = setdiff((1:n),badI); % bad estimation/data indexes
good_samp_idx = goodI(randperm(length(goodI),num2disp));
good_samp = est(good_samp_idx);
if ~isempty(badI)
%     bad_samp_idx = badI(randperm(length(badI),2));
    bad_samp = est(badI(1:num2disp));
    flag = 1;
end
%% Display numerous examples

% Display good samples of the estimated acorr_CGP and the real data
for m = 1:num2disp
        plotFit( X, good_samp(m).Results.Data , good_samp(m).Results.estimation, 'Simulated kernel fit test' )
end

% Display bad samples of the estimated acorr_CGP and the real data
if flag
    for m = 1:num2disp
            figure;
            plot(X,bad_samp(m).Results.estimation,' +r','LineWidth',0.4);hold on;
            plot(X, bad_samp(m).Results.Data,'-b');
            hold off;
            figTitle = ['diff in parameters = ',num2str(bad_samp(m).PEresults.paramData - bad_samp(m).PEresults.paramEst)];
            plotFit( X, bad_samp(m).Results.Data , bad_samp(m).Results.estimation, figTitle)
            title();
    end;
end;
%% Analysis of entire Data measurments

Rsquare_bar = mean(Rsquare);
Rsquare_std = std(Rsquare);
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