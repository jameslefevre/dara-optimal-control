% Blocks of code in this script control runs of the optimal control algorithm. 
% Parameters are customised, the initial state determined, and
% the bang-bang or continuous optimal control function is called.
% The parameter saveString controls the folder and filenames where results
% are stored. 
% Loops are not included here but are advised when running multiple parameter combinations.


% example: continuous and bang-bang convergence for full default model

clear
params = getparams();
params("omega") = 0.9;
params("saveString") = "v0_cont_w9";
APN_vals = steadystates_dara_immune(params,0,1e-24); 
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 1000;
params("MinIters") = 200;
params("iterationsPlot") = [10,50,100,150,200,400,600,800,1000];

[t_y,y,U,iterations,converged,ConvergenceStats] = contcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})


clear
params = getparams();
params("omega") = 0.02;
params("saveString") = "v0_w02";
APN_vals = steadystates_dara_immune(params,0,1e-24); 
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 100000;
params("MinIters") = 10000;
params("iterationsPlot") = [10,100,1000,5000,10000,20000,30000,50000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})

% example: continuous and bang-bang convergence for null-N model 
% (control model without drug escape and off-target effect)

clear
params = getparams();
params("pp") = 0.27; %Slightly reduced proliferation of P
params("mp") = 0.05; %Slightly increased exit of P
params("dp") = 0.0; %Disable loss of CD38 expression in P
params("pn") = 0.0; %Disable proliferation of N
params("mn") = 1.0; 
params("dn") = 1.0; %Rapid death or gain of CD38 expression in any N present
params("mau") = 0.0; %0.5; %Additional exit of A per unit of control; this disables off-target effect
params("dpu") = 0.0; %1.0; %Additional loss of CD38 expression in P per unit of control; this disables loss of CD38 expression due to drug

params("omega") = 0.9;
params("saveString") = "nullN_cont_w9";
APN_vals = steadystates_dara_immune(params,0,1e-24); 
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 1000;
params("MinIters") = 200;
params("iterationsPlot") = [10,50,100,150,200,400,600,800,1000];

[t_y,y,U,iterations,converged,ConvergenceStats] = contcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})


clear
params = getparams();
params("pp") = 0.27; %Slightly reduced proliferation of P
params("mp") = 0.05; %Slightly increased exit of P
params("dp") = 0.0; %Disable loss of CD38 expression in P
params("pn") = 0.0; %Disable proliferation of N
params("mn") = 1.0; 
params("dn") = 1.0; %Rapid death or gain of CD38 expression in any N present
params("mau") = 0.0; %0.5; %Additional exit of A per unit of control; this disables off-target effect
params("dpu") = 0.0; %1.0; %Additional loss of CD38 expression in P per unit of control; this disables loss of CD38 expression due to drug

params("omega") = 0.02;
params("saveString") = "nullN_w02";
APN_vals = steadystates_dara_immune(params,0,1e-24); 
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 100000;
params("MinIters") = 10000;
params("iterationsPlot") = [10,100,1000,5000,10000,20000,30000,50000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})



% example: simple bang-bang convergence case

clear
params = getparams();
params("mau") = 0.5;
params("dp") = 0.0003;
params("dn") = 0.003;
params("dpu") = 0.2;
% params("alpha_") = 0;
params("omega") = 0.02;
params("saveString") = "grid_2_2_3_w02";
APN_vals = steadystates_dara_immune(params,0,1e-24); % calculate the steady state solutions
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 100000;
params("MinIters") = 10000;
params("iterationsPlot") = [10,100,1000,5000,10000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})

% example: intermediate bang-bang convergence case (note reduced omega)

clear
params = getparams();
params("mau") = 0.5;
params("dp") = 0.0003;
params("dn") = 0.003;
params("dpu") = 2.0;
params("alpha_") = 0;
params("omega") = 0.01;
params("saveString") = "grid_2_20_3_0_w01";
APN_vals = steadystates_dara_immune(params,0,1e-24); % calculate the steady state solutions
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 100000;
params("MinIters") = 20000;
params("iterationsPlot") = [10,100,1000,5000,10000,15000,20000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})

% example: one of the most difficult bang-bang convergence cases 

clear
params = getparams();
params("mau") = 0.5;
params("dp") = 0.0003;
params("dn") = 0.003;
params("dpu") = 0.2;
params("alpha_") = 0;
params("omega") = 0.02;
params("saveString") = "grid_2_2_3_0_w02";
APN_vals = steadystates_dara_immune(params,0,1e-24); % calculate the steady state solutions
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

params("MaxIters") = 100000;
params("MinIters") = 0;
params("iterationsPlot") = [10,100,1000,5000,10000,20000,30000,40000,50000,60000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})



% example: continuing an optimal control calculation

% Continue the bang-bang convergence calculation for default model. 
% This code requires the original calculation to be run so that the result
% at t=10000 is available to be reloaded and continued.
% This case is fully converged, but code is designed for complex cases
% in which the calculation has been interrupted prematurely.

clear
params = getparams();
params("omega") = 0.02;
params("saveString") = "v0_w02";
APN_vals = steadystates_dara_immune(params,0,1e-24); 
params('A_init') = APN_vals{end,'A'};
params('P_init') = APN_vals{end,'P'};
params('N_init') = APN_vals{end,'N'};

open('v0_w02\v0_w02_10000.fig'); % data at checkpoints saved in fig files
ydata = get(get(gca,'Children'), 'YData'); % order reversed vs legend in the fig
params('iteration_init') = 10000; % restart t where we left off
params('U_init') = ydata{1}; % The info we need to continue is the u solution.
close(gcf);

params("MaxIters") = 100000;
params("MinIters") = 11000;
params("iterationsPlot") = [10001,11000];

[t_y,y,U,iterations,converged,ConvergenceStats] = bbcontrol_dara_immune(params);
disp(ConvergenceStats{end,'Cost'})






