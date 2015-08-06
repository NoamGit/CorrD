function [x,fval,exitflag,output] = annealData2Acorr(x0,lb,ub, costfunction)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = saoptimset;
%% Modify options setting
options = saoptimset(options,'TemperatureFcn', @temperatureboltz);
options = saoptimset(options,'AcceptanceFcn', @custom_acceptancesa);
options = saoptimset(options,'Display', 'off');
options = saoptimset(options,'HybridFcn', {  @patternsearch [] });
[x,fval,exitflag,output] = ...
simulannealbnd(@(p)costfunction(p),x0,lb,ub,options);