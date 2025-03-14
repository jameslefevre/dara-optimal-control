function [params] = getparams()
% Create a map of parameters with default values

params = containers.Map;
params("Tfinal") = 200;  %Specified final time
params("dt") = 0.001; %Time-step
params("omega") = 0.1; %Convergence parameter, 0.9 gives good result for cont control
params("RelTol") = 1e-3; %Desired relative tolerance for convergence
params("RelTol1000") = 1e-2; %Desired relative tolerance over 1000 time steps for convergence
params("MaxIters") = 250; %Number of iterations to perform before giving up if convergence is not reached
params("MinIters") = 0; % force at least this many iterations regardless of convergence (but MaxIters overrides)
params("iterationsPlot") = []; % save plots at these times; allows visual tracking on convergence, and simulation may also be restarted from a .fig file
params("saveString") = join("unnamed"); % name of folder to save results; filenames also contain this string

%Model parameters
params("ba") = 0.1008; %Influx of A % had 0.72
params("pa") = 0.43; %Proliferation of A
params("ma") = 0.44; %Exit of A
params("pp") = 0.28; %Proliferation of P
params("mp") = 0.048; %Exit of P
params("dp") = 0.003; %Loss of CD38 expression in P
params("pn") = 0.15; %Proliferation of N
params("mn") = 0.06; %Exit of N
params("dn") = 0.03; %Gain of CD38 expression in N

params("alpha_") = 0.015; % immune response rate ; values from Sharp20019
params("gamma_") = 0.1; % half saturation constant of immune response % table 1 says 0.01, fig 4 says 0.1

params("mau") = 0.1; %0.5; %Additional exit of A per unit of control
params("mpu") = 1.0; %Additional exit of P per unit of control
params("dpu") = 0.2; %1.0; %Additional loss of CD38 expression in P per unit of control

%Pay-off weightings, default values
params("a1") = 1; %Weighting on negative impact of chemotherapy control
params("a2") = 1; %Weighting on negative impact of leukaemia

%Control bounds determining minimum (lower) and maximum (upper) admissable bounds of
%the chemotherapy bang-bang control.
params("Ulower") = 0; 
params("Uupper") = 1.0; 

% set default initial values to 0 - these are typically set using steady state
% solution, depends on other parameters
params("A_init") = 0; 
params("P_init") = 0; 
params("N_init") = 0; 

params("U_init") = 0; % can use this to set a non-trivial starting point for the control U
% if length 1, the value will be applied at all time steps; otherwise require length Nt 
params("iteration_init") = 0; % use with U_init to continue an optimal control convergence process;
% first iteration of new run will be numbered iteration_init+1

