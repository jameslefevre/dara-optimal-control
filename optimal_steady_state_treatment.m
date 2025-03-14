function [cost,U,A,P,N]  = optimal_steady_state_treatment(params,CostFn, controlGrid, tol)
% optimal_steady_state_treatment looks for the value of u which minimises the
% provided cost function (CostFn), assuming a steady state of the system
% parameterised by params
% Uses steady_state_cost function to do the main calculation,
% and system function fminbnd to get min.
% Starts with grid search over provided controlGrid to minimise chance of
% getting local min; optimal u is assumed to lie within the range of
% controlGrid, but it is easy to verify the solution if the cost is less
% than the control-only cost of the maximum u value in controlGrid 




fn = @(u) [steady_state_cost(params,CostFn,u, tol)];

% tol=1e-24;
%fn = @(u) [6*(u-0.5)^4 - (u-0.5)^2]; % testing

sz = size(controlGrid);
gridCosts = zeros(sz);
n = sz(1,2);
for i = 1:n
    gridCosts(1,i) = fn(controlGrid(i));
end

gridMins = [];
if gridCosts(1,1) < gridCosts(1,2)
    gridMins(end+1) = 1;
end
for i = 2:(n-1)
    if (gridCosts(1,i) <= gridCosts(1,i-1)) && (gridCosts(1,i) <= gridCosts(1,i+1)) 
        gridMins(end+1) = i;
    end
end
if gridCosts(1,n) < gridCosts(1,n-1)
    gridMins(end+1) = n;
end

uMins = [];
cMins = [];
for i = gridMins
    a = controlGrid(max(i-1,1));
    b = controlGrid(min(i+1,n));
    [uMins(end+1),cMins(end+1)] = fminbnd(fn,a,b);
end
[~,ind] = min(cMins);
controlLevel=uMins(ind);
[cost,U,A,P,N] = steady_state_cost(params,CostFn,controlLevel, tol);
end