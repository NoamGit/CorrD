function [ ] = plotFit( x, y , yfit, str_title )
% Create a scatter plot of the original x and y data
figure
s = 3;
scatter(x, y, s, 'MarkerEdgeColor','c')

% Plot yfit
line(x, yfit, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 1)

% Add a legend and axis labels
legend('Data', 'Fit', 'Location', 'SouthEast')
xlabel('Time [s]')
ylabel('Correlation')
title(str_title);
axis tight;

% % Build a string that contains the Latex expression
eqtext = '\itf(t)={\itexp({{-(t-\mu)^2}/{\sigma^2}})}';
eqtext = [eqtext '\itsin({2\pif(t-\phi)})'];
% % Add the string containing the Latex expression to the plot
text(x(ceil(numel(x)/8)), max(y)/(-2), eqtext, 'FontSize', 12, 'Color', 'k')
end

