function [] = plotEstimation( obj, data )
%plotEstimation( obj, data ) Plots Parameter estimation given the data
    
    % Define Variables
    %     X = data.linKernel.timelag;
    X = obj.acorr.lags;
    Ynfilt = obj.acorr.CGP_Corr;
    Y_forw = data.acorr.CGP_Corr;
    Y_est = obj.PEresults.OptimResults.estimation;
    goodness = obj.PEresults.GoodnessStatistics;
    
    hAll = figure;
    set(hAll, 'Position', [100, 150, 1100, 650]);
    hold on
    hRawData_back = line(X, normax(Ynfilt));
    hData_forw = line(X, normax(Y_forw));
    hFit = line(X  , normax(Y_est));

    % Adjust line properties (functional)
    set(hFit, 'Color', [0 0 .5])
    set(hRawData_back, 'LineStyle', 'none', 'Marker', '.')
    set(hData_forw, 'LineStyle', '-.', 'Color', 'r') % VV

    % Adjust line properties (aesthetics)
    set(hFit, 'LineWidth', 2) % VV
    set(hRawData_back, 'Marker', 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1]) % VV
    set(hData_forw, 'LineWidth', 1.5)

    % Add labels
    hTitle = title('Acorr comparsion estimated vs. real data');
    hXLabel = xlabel('Time [s]');
    hYLabel = ylabel('normalized Correlation');

    % Add legend
    hLegend = legend([hFit, hRawData_back, hData_forw], ...
        'Estimation {e^{(-\itx^2)/2\it\sigma^2}}cos({2\it\pif\cdot x}) - e^{(-\it x^2 - (2\it\pi f)^2\it\sigma^4)/{2\it\sigma^2}}',...
        'CGP acorr back process', ... 
        'CGP acorr forward process', ...
        'Location', 'NorthWest');

    % Adjust font
    set(gca, 'FontName', 'Helvetica')
    set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
    set([hLegend, gca], 'FontSize', 10)
    set([hXLabel, hYLabel], 'FontSize', 14)
    set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
    
    % Create textbox
    annotation(hAll,'textbox',[0.65625 0.776483050847458 0.110416666666667 0.0656779661016949],...
        'String',{['\itR^2 =  ',num2str(goodness.rsquare)] ,['\itSSE = ',num2str(goodness.sse)]},'FitBoxToText','on',...
    'EdgeColor','none','FontSize',16);

    % Adjust axes properties
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3]);
    axis tight;
    hold off;
end

