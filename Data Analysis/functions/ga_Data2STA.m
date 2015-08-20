function [x,fval,exitflag,output,population,score] = ga_Data2STA(nvars,lb,ub,costfunction)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', [-1;1]); % like in MATLAB 2013
options = gaoptimset(options,'PopulationSize', 200); 
options = gaoptimset(options,'TolCon', 1e-6); % like in MATLAB 2013
options = gaoptimset(options,'HybridFcn', {  @patternsearch [] });
options = gaoptimset(options,'Display', 'off');
[x,fval,exitflag,output,population,score] = ...
ga(@(p)costfunction(p),nvars,[],[],[],[],lb,ub,[],[],options);
