%% Report for estSTA_Sim_3_1
% run first estSTA_Data_3_1

load('D:\# Projects (Noam)\# CorrD\code\trunk\Workspace\WS_estSTA_3_1_X20.mat');

Rsquare = cell2mat({goodness.rsquare});
estSTA = optimResults.estimationSTA;
n = numel(realSTA);

flag = 0;
num2disp = 4;
badI = find(Rsquare < 0.85); % good estimation/data indexes
goodI = setdiff((1:n),badI); % bad estimation/data indexes
good_samp_idx = goodI(randperm(length(goodI),num2disp));
good_samp = estSTA(good_samp_idx, :);
if ~isempty(badI)
%     bad_samp_idx = badI(randperm(length(badI),2));
    bad_samp = estSTA(badI(1:num2disp), :);
    flag = 1;
end
%% Display numerous examples

% Display good samples of the estimated acorr_CGP and the real data
for m = 1:num2disp
        plotFit( time, normax(realSTA{good_samp_idx(m)}) , good_samp(m,:), 'Real Kernel fit test good samples')
end

% Display bad samples of the estimated acorr_CGP and the real data
if flag
    for m = 1:num2disp
%             figure;
%             plot(X,bad_samp(m).Results.estimation,' +r','LineWidth',0.4);hold on;
%             plot(X, bad_samp(m).Results.Data,'-b');
%             hold off;
%             figTitle = ['diff in parameters = ',num2str(bad_samp(m).PEresults.paramData - bad_samp(m).PEresults.paramEst)];
            plotFit( time, normax(realSTA{badI(m)}) , bad_samp(m,:), 'Real Kernel fit test bad samples');           
    end
end
%% Analysis of entire Data measurments

Rsquare = cell2mat({goodness.rsquare});
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