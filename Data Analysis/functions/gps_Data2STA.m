function [x,fval,exitflag,output] = gps_Data2STA(x0,lb,ub,costfunction)
%% This is an auto generated MATLAB file from Optimization Tool.
% check if is same as in MATLAB 2013Bb

%% Start with the default options
options = psoptimset;
%% Modify options setting
options = psoptimset(options,'Display', 'off');
[x,fval,exitflag,output] = ...
patternsearch(@(p)costfunction(p),x0,[],[],[],[],lb,ub,[],options);
