%% grid search for equilibria, default model parameters
% numerical check that there are two stable equilibria

clear
params = getparams();
%controlGrid = 0:0.01:1;
controlGrid = 0:0.1:1;

APN_vals = table('Size',[length(controlGrid)^3,6],'VariableTypes',["double","double","double","double","double","double"] , 'VariableNames',["A0","P0","N0","A","P","N"] );
i = 0;
for A0 = controlGrid
    for P0 = controlGrid
        for N0 = controlGrid
            tic
            [A,P,N] = steadystate_from_specified_start_dara_immune(params,A0,P0,N0,0,1e-10);
            i=i+1;
            APN_vals{i,"A0"} = A0;
            APN_vals{i,"P0"} = P0;
            APN_vals{i,"N0"} = N0;
            APN_vals{i,"A"} = A;
            APN_vals{i,"P"} = P;
            APN_vals{i,"N"} = N;
            [A0,P0,N0,A,P,N]
            toc
        end
    end
end


[min(APN_vals.A),max(APN_vals.A)] % 0.3180, 0.4727
[min(APN_vals.P),max(APN_vals.P)] % 0, 0.3813
[min(APN_vals.N),max(APN_vals.N)] % 0, 0.0148

delt = 0.005; %1e-3;
[min(APN_vals.A(APN_vals.P>delt)),max(APN_vals.A(APN_vals.P>delt))] % 0.3180
[min(APN_vals.P(APN_vals.P>delt)),max(APN_vals.P(APN_vals.P>delt))] % 0.3813
[min(APN_vals.N(APN_vals.P>delt)),max(APN_vals.N(APN_vals.P>delt))] % 0.0148

[min(APN_vals.A(APN_vals.P<=delt)),max(APN_vals.A(APN_vals.P<=delt))] % 0.4727
[min(APN_vals.P(APN_vals.P<=delt)),max(APN_vals.P(APN_vals.P<=delt))] % 0
[min(APN_vals.N(APN_vals.P<=delt)),max(APN_vals.N(APN_vals.P<=delt))] % 0

mean(APN_vals.P>delt) % 0.9812 
sum(APN_vals.P<=delt) % 25

% this grid only gives the N=P=0 solution when P0=0; there are more cases
% where P \in (0,0.1) though, e.g.
[A,P,N] = steadystate_from_specified_start_dara_immune(params,0.8,0.05,0.05,0,1e-10);
[A,P,N] 

%% infinite horizon optimisation
% two zero and two non-zero examples

% default params, linear cost 
params = getparams();
controlGrid = 0:0.01:1;
CostFn = @(U,PN) [U+PN];
[cost,U,A,P,N]  = optimal_steady_state_treatment(params,CostFn, controlGrid, 1e-24);
"linear, default"
[cost,U,A,P,N]

% default params, quadratic cost 
params = getparams();
controlGrid = 0:0.01:1;
CostFn = @(U,PN) [U*U+PN*PN];
[cost,U,A,P,N]  = optimal_steady_state_treatment(params,CostFn, controlGrid, 1e-24);
"quadratic, default"
[cost,U,A,P,N]

% no immune, linear cost 
params = getparams();
params("alpha_")=0;
controlGrid = 0:0.01:1;
CostFn = @(U,PN) [U+PN];
[cost,U,A,P,N]  = optimal_steady_state_treatment(params,CostFn, controlGrid, 1e-24);
"linear, no immune"
[cost,U,A,P,N]

% no immune, quadratic cost 
params = getparams();
params("alpha_")=0;
controlGrid = 0:0.01:1;
CostFn = @(U,PN) [U*U+PN*PN];
[cost,U,A,P,N]  = optimal_steady_state_treatment(params,CostFn, controlGrid, 1e-24);
"quadratic, no immune"
[cost,U,A,P,N]

