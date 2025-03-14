function [t_y,y,U,iterations,converged,ConvergenceStats] = generalcontrol_dara_immune(params,kappa)
%generalcontrol_dara_immune Calculates the optimal control for the
% Darra/MM model with immune response, and cost function
% u+P+N+kappa(u^2+(P+N)^2)
% The control u is constrained by the bounds Ulower, Uupper given in params

%   Uses iterative back/forward sweep method with state/costate solving
%   2-point BVP. 

% export parameters from params map to local namespace
for key = params.keys()
    eval(append(key{1}," = params('",key{1},"');"));
end
Nt = Tfinal/dt+1;
t_y = linspace(0,Tfinal,Nt);

% folder for results 
if ~isfolder(saveString)
    mkdir(saveString);
end

% initialise
% have not implemented U_init stuff here; used in bbcontrol_dara_immune to allow resumptions
y(:,1) = [A_init,P_init,N_init];

U = zeros(1,Nt);
uold =  U;
iterations = 0; %initialise iteration count
converged = false;
RelTolTest = 0;
%Initialise vectors to monitor convergence
RelativeTolerance = []; 
SumU = [];
SumPN = [];
SumU2 = [];
SumPN2 = [];

%State equations
State = @(t,y,U) [(ba+pa*y(1)*(1-y(1)-y(2)-y(3))-ma*y(1)-mau*U*y(1)); 
    (pp*y(2)*(1-y(1)-y(2)-y(3)) -dp*y(2) +dn*y(3) -dpu*U*y(2) -mp*y(2)-mpu*U*y(2) - alpha_*y(2)/(gamma_ + y(2)+y(3)) ); 
    (pn*y(3)*(1-y(1)-y(2)-y(3)) +dp*y(2) -dn*y(3) +dpu*U*y(2) -mn*y(3)) - alpha_*y(3)/(gamma_ + y(2)+y(3)) ];

%Costate equations ; a2 term in part 2 and 3 is cost function contribution
Costate = @(t,Lambda,y,U) [Lambda(1)*(-pa*(1-2*y(1)-y(2)-y(3))+ma+mau*U) + Lambda(2)*pp*y(2) + Lambda(3)*pn*y(3);
	-1-2*kappa*(y(2)+y(3))*a2  + Lambda(1)*pa*y(1) + Lambda(3)*(pn*y(3)-dp-dpu*U-alpha_*y(3)/(gamma_+y(2)+y(3))^2 ) + Lambda(2)*(-pp*(1-y(1)-2*y(2)-y(3))+dp+dpu*U+mp+mpu*U + alpha_*(gamma_+y(3))/(gamma_+y(2)+y(3))^2 );
	-1-2*kappa*(y(2)+y(3))*a2 + Lambda(1)*pa*y(1) + Lambda(2)*(pp*y(2)-dn-alpha_*y(2)/(gamma_+y(2)+y(3))^2) + Lambda(3)*(-pn*(1-y(1)-y(2)-2*y(3))+dn+mn+ alpha_*(gamma_+y(2))/(gamma_+y(2)+y(3))^2)]; 

% Forward-backward sweep
while(iterations<MaxIters)
      
i = 0; %Initialise loop variable
t = 0; %Initialise time for forward sweep

%Forward sweep using fourth-order Runge-Kutta scheme
% disp(Nt);disp(size(t));disp(size(y));disp(size(U))
while i < Nt-1
    t = t + dt;
    i=i+1;
    k1 = State(t,y(:,i),U(i));
    k2 = State(t,y(:,i)+dt*k1/2,0.5*(U(i)+U(i+1)));
    k3 = State(t,y(:,i)+dt*k2/2,0.5*(U(i)+U(i+1)));
    k4 = State(t,y(:,i)+dt*k3,U(i+1));
    y(:,i+1) = y(:,i) + (dt/6)*(k1+2*k2+2*k3+k4);
end

t = Tfinal; %Initialise time for backward sweep
Lambda = zeros(3,length(y)); %Initialise Lambda
Lambda(:,end) = [0,0,0]; %Apply transversality conditions to obtain final time condition on Lambda (costate)
j = Nt; %Initialise loop variable

%Backward sweep using fourth-order Runge-Kutta scheme 
while j > 1 
    t = t - dt;
    j=j-1;
    k1 = Costate(t,Lambda(:,j+1),y(:,j+1),U(j+1));
    k2 = Costate(t,Lambda(:,j+1)-dt*k1/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)));
    k3 = Costate(t,Lambda(:,j+1)-dt*k2/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)));
    k4 = Costate(t,Lambda(:,j+1)-dt*k3,y(:,j),U(j));
    Lambda(:,j) = Lambda(:,j+1) - (dt/6)*(k1+2*k2+2*k3+k4);
end

% record all stats / print results here, before updating U
% except for RelativeTolerance, which is recorded to the stats after U update

sumU = sum(U)*dt;
sumPN = sum(y(2,:)+y(3,:))*dt;
sumU2 = sum(U.^2)*dt;
sumPN2 = sum((y(2,:)+y(3,:)).^2)*dt;
SumU = [SumU ; sumU];
SumPN = [SumPN ; sumPN];
SumU2 = [SumU2 ; sumU2];
SumPN2 = [SumPN2 ; sumPN2];


%Update and display iteration count
iterations = iterations+1;
fprintf('Rel tol: %d  ;  ',RelTolTest) 
fprintf('Sum U^2,(P+N)^2,cost: %d,%d,%d;  ',sumU2,sumPN2,a1*sumU2+a2*sumPN2)
fprintf('Iters: %d  \n',iterations) 

if ismember(iterations,iterationsPlot)
     fig = plot_optimal_control_and_variables_over_time(containers.Map({'t_y','y','U','plotPNsum','saveName','includeLegend'},...
           {t_y,y,U,false, append(saveString,", ",num2str(iterations)),true}));
     saveas(fig,append(saveString,"/",saveString,'_',num2str(iterations),'.fig') );
     saveas(fig,append(saveString,"/",saveString,'_',num2str(iterations),'.png') );
     close(fig);
end

uold =  U; %Store control from previous iteration

Uupdate = ( Lambda(1,:).*mau.*y(1,:) +Lambda(2,:).*dpu.*y(2,:) +Lambda(2,:).*mpu.*y(2,:) - Lambda(3,:).*dpu.*y(2,:) - 1) ./ (2*a1*kappa);
U = omega*U + (1-omega)*Uupdate;
U = min(max(U,Ulower),Uupper);

%Check for convergence
RelTolTest = RelTol*norm(U) - norm(U-uold);
RelativeTolerance = [RelativeTolerance ; RelTolTest];

if (iterations>MinIters && RelTolTest > 0) 
    fprintf('Specified relative tolerance of %g has been met \n\r',RelTol)
    converged = true;
    break
end  

end

% final sequence: save final plot, stats, convergence plot 

fig = plot_optimal_control_and_variables_over_time(containers.Map({'t_y','y','U','plotPNsum','saveName','includeLegend'},...
           {t_y,y,U,false, append(saveString,", ",num2str(iterations)),true}));
saveas(fig,append(saveString,"/",saveString,'_',num2str(iterations),'_fin.fig') );
close(fig);

Cost = a1 .* SumU + a2 .* SumPN + kappa .* a1 .* SumU2 + kappa .* a2 .* SumPN2;

ConvergenceStats = table(SumU,SumPN,SumU2,SumPN2,Cost,RelativeTolerance);
writetable(ConvergenceStats,append(saveString,"/",saveString,"_stats.csv"));

colours = [ 
    0/255  114/255  189/255
    222/255  125/255  0/255 
]; % Define colours for convergence plot

iters = 1:length(Cost);
fig = figure('Name',append(saveString," convergence") );
ax = axes('Parent', fig);
set(ax, 'ColorOrder', colours);
hold on
box on
grid on
line1 = plot(iters,Cost,'LineWidth',2);
line2 = plot(iters,RelativeTolerance,'LineWidth',2);
legend([line1,line2],{'Cost','RelativeTolerance'},'Location','northeast');
xlabel('Iterations','fontsize',18);
%axis([0,t_y(end),0,1]) % axis([0,Tfinal,0,1])
%xt = get(ax, 'XTick');
set(ax, 'FontSize', 18)
saveas(fig,append(saveString,"/",saveString,'_convergence.fig') );
saveas(fig,append(saveString,"/",saveString,'_convergence.png') );
close(fig);

