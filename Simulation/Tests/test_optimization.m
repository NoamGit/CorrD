%% test the optimization script
% checks if the optimization is done right

[ est, proc, filt_coeff ] = test_generateDummyProc( 3 );
%% Single process Isolation
idx = 3;
X = est{idx}.acorr.lags;
Ynfilt = est{idx}.acorr.CGP_Corr;
Y_back = est{idx}.acorr.filt_acorr;
Y_forw = proc{idx}.acorr.CGP_Corr;
parameters = proc{idx}.linKernel.parameters([1,4]); % [sig1, ~, ~, f]
N = length(est{idx}.acorr.lags);

% plot compare Y back, Y forw,Y not filt
plot(X,[Y_back Y_forw Ynfilt]);
legend('Y back', 'Y forw','Y not filt')
    %% Plot - computing average power spectrum of signal to find dominant freq ( = f)
    fest = fft(est{idx}.acorr.CGP_Corr); % single
    Pyy_est = fest.*conj(fest)/N;
    f = est{idx}.Fs/N*(0:N/2);
    figure;
    plot(f,Pyy_est(1:N/2+1))
    title('Power spectral density')
    xlabel('Frequency (Hz)')

%% Optimization
    %% fit 
    
    Lower = [eps 0 0]; % [ N f sig ]
    Upper = [Inf 8 1];
    modelEq = 'N * (exp((-x.^2)/(2*s^2)).*cos((2*pi*f).*x))-exp((-x.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0)';
    opt_fit{idx} = paramEstim(est{idx}, proc{idx}, X, modelEq , Lower, Upper, filt_coeff{idx}  );
    opt_fit_re = opt_fit{idx}.PEresults.fitResults;
    opt_fit_stat = opt_fit{idx}.PEresults.fitsStatistics;
    %% Optimization

    Lower = [eps 0 0]; % [ N f sig ]
    Upper = [Inf 10 10];
    fun = @(N,f,s) N * (exp((-X.^2)/(2*s^2)).*cos((2*pi*f).*X))-exp((-X.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0);
    costfun = @(p) sum( (Y_back' - fun(p(1),p(2),p(3))).^2 ); % Non Linear Least Square 
    x0 = [1 1 1];
        %% Plots
        % parametric space visulization - determine Nt
        Nt = 0.179;
        ezfun = @(f,s) sum( (Y_back' - Nt * (exp((-X.^2)/(2*s^2)).*cos((2*pi*f).*X))-exp((-X.^2-(2*pi*f)^2*s^4)/(2*s^2)).*cos(2*(2*pi*f)*0)).^2 );
        figure;
        ezsurf(ezfun, [0 10 0 10])

        % check model with optial parameters 
        figure;
        plot(X, Y_back,X , fun(1.2, parameters(2), parameters(1)), X, 5*Y_forw );
        legend('Y back', 'Y optimal param', 'Y forw')
    %% SA - slow but works for most cases
    [x,fval,exitflag,output] = annealData2Acorr(x0, Lower,Upper,costfun);
    Y_est = fun( x(1), x(2), x(3));
    numcoeff = 3;
    dfe = N - numcoeff; % degrees of freedom - number of observations - num of coefficients
    res = Y_back' - Y_est;
    goodness = iGoodnessStructure(Y_est,[],res,dfe,N);
    
    %% Plot
    % Create basic plot
    
    hAll = figure;
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
    
    % Add text
%     hText = text(10, 800, ...
%         sprintf('{\itR^2 = %d   \itSSE = %d}', goodness.rsquare, goodness.sse));
    % Create textbox
    annotation(hAll,'textbox',[0.65625 0.776483050847458 0.110416666666667 0.0656779661016949],...
        'String',{['\itR^2 =  ',num2str(goodness.rsquare)] ,['\itSSE = ',num2str(goodness.sse)]},'FitBoxToText','off',...
    'EdgeColor','none','FontSize',16,'FontName', 'Helvetica');

    % Adjust axes properties
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3]);