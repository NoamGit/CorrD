function acceptpoint = custom_acceptancesa(optimValues,newx,newfval)    
%ACCEPTANCESA Acceptance function for simulated annealing solver
%   ACCEPTPOINT = ACCEPTANCESA(optimValues,newX,newfval) uses the
%   change in function values between the current point and new point to
%   determine whether the new point is accepted or not.
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter
%
%   NEWX: new point 
%
%   NEWFVAL: function value at NEWX
%
%   Example:
%    Create an options structure using ACCEPTANCESA as the annealing
%    function
%    options = saoptimset('AcceptanceFcn',@acceptancesa);

%   Copyright 2006-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:21:27 $

% new parameters
G = 2000; % After G iteration the probability to accept a step is P2 (before it is P1)
P1 = 0.35;
P2 = 0.1;

delE = newfval - optimValues.fval;

% If the new point is better accept it
if delE < 0
    acceptpoint = true;
% Otherwise, accept it randomly based on a Boltzmann probability density
else
    h = 1/(1+exp(delE/max(optimValues.temperature)));
    if optimValues.iteration < G  
        if h > P1
            acceptpoint = true;
        else
            acceptpoint = false;
        end
    else
         if h > P2
            acceptpoint = true;
        else
            acceptpoint = false;
         end
    end
end