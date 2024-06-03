% examples of simple simulations (no optimal control)


% first example: default and null-N models with and without immune response
% no control, t=50

clear
colours = [ 
    30/255  136/255  229/255 % blue; cost (and A)
    216/255  27/255  96/255 % red; P
    255/255  193/255  7/255 % yellow; N  
];

for noImmune = [false,true]
    for nullN = [false,true]
        params = getparams();     
        params("Tfinal") = 50;
        if noImmune
            params("alpha_") = 0.0;
            nm = "simulation_noImmune";
        else
            nm = "simulation_immune";
        end
        if nullN
            nm = join([nm,"_nullN"],"");
            params("pp") = 0.27; %Slightly reduced proliferation of P
            params("mp") = 0.05; %Slightly increased exit of P
            params("dp") = 0.0; %Disable loss of CD38 expression in P
            params("pn") = 0.0; %Disable proliferation of N
            params("mn") = 1.0; 
            params("dn") = 1.0; %Rapid death or gain of CD38 expression in any N present
            params("mau") = 0.0; %0.5; %Additional exit of A per unit of control; this disables off-target effect
            params("dpu") = 0.0; %1.0; %Additional loss of CD38 expression in P per unit of control; this disables loss of CD38 expression due to drug
        end
        for key = params.keys()
            eval(append(key{1}," = params('",key{1},"');"));
        end
        Nt = Tfinal/dt+1;
        t_y = linspace(0,Tfinal,Nt);

        y = zeros(3,Nt);
        y(:,1) = [0,0.1,0]; % 2a,4a - A=0,L=0.1
        U = zeros(1,Nt);

        %State equations
        State = @(t,y,U) [(ba+pa*y(1)*(1-y(1)-y(2)-y(3))-ma*y(1)-mau*U*y(1)); 
            (pp*y(2)*(1-y(1)-y(2)-y(3)) -dp*y(2) +dn*y(3) -dpu*U*y(2) -mp*y(2)-mpu*U*y(2) - alpha_*y(2)/(gamma_ + y(2)+y(3)) ); 
            (pn*y(3)*(1-y(1)-y(2)-y(3)) +dp*y(2) -dn*y(3) +dpu*U*y(2) -mn*y(3)) - alpha_*y(3)/(gamma_ + y(2)+y(3)) ];
        
        %Forward sweep using fourth-order Runge-Kutta scheme
        i = 0; %Initialise loop variable
        t = 0; %Initialise time for forward sweep
        while i < Nt-1
            t = t + dt;
            i=i+1;
            k1 = State(t,y(:,i),U(i));
            k2 = State(t,y(:,i)+dt*k1/2,0.5*(U(i)+U(i+1)));
            k3 = State(t,y(:,i)+dt*k2/2,0.5*(U(i)+U(i+1)));
            k4 = State(t,y(:,i)+dt*k3,U(i+1));
            y(:,i+1) = y(:,i) + (dt/6)*(k1+2*k2+2*k3+k4);
        end

        figure
        set(gca, 'ColorOrder', colours);
        hold on
        box on
        line1 = plot(t_y,y(1,:),'LineWidth',2);
        line2 = plot(t_y,y(2,:),'LineWidth',2);
        line3 = plot(t_y,y(3,:),'LineWidth',2);
        hL = legend([line1,line2,line3],{'A','P','N'},'Location','northeast');
        ylabel('State','fontsize',18);
        xlabel('Time','fontsize',18);
        axis([0,Tfinal,0,1])
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 18)
        set(gcf, 'Position',  [100, 100, 400, 480]);

        saveas(gcf,append(nm,'.png') );
        close(gcf);

    end
end



% second example: default model and null-N model with two predefined
% controls; t=50

clear
colours = [ 
    30/255  136/255  229/255 % blue; cost (and A)
    216/255  27/255  96/255 % red; P
    255/255  193/255  7/255 % yellow; N  
];

for controlCase = [0,1]
    for nullN = [false,true]
        params = getparams();  
        params("Tfinal") = 50;
        if controlCase==0
            nm = "simulation_control1";
        else
            nm = "simulation_control2";
        end
        if nullN
            nm = join([nm,"_nullN"],"");
            params("pp") = 0.27; %Slightly reduced proliferation of P
            params("mp") = 0.05; %Slightly increased exit of P
            params("dp") = 0.0; %Disable loss of CD38 expression in P
            params("pn") = 0.0; %Disable proliferation of N
            params("mn") = 1.0; 
            params("dn") = 1.0; %Rapid death or gain of CD38 expression in any N present
            params("mau") = 0.0; %0.5; %Additional exit of A per unit of control; this disables off-target effect
            params("dpu") = 0.0; %1.0; %Additional loss of CD38 expression in P per unit of control; this disables loss of CD38 expression due to drug
        end
        for key = params.keys()
            eval(append(key{1}," = params('",key{1},"');"));
        end
        Nt = Tfinal/dt+1;
        t_y = linspace(0,Tfinal,Nt);

        APN_vals = steadystates_dara_immune(params,0,1e-14); % calculate the steady state solutions
        y = zeros(3,Nt);
        y(:,1) = [APN_vals{end,'A'},APN_vals{end,'P'},APN_vals{end,'N'}];
        U = zeros(1,Nt);
        if controlCase==0
            U = zeros(1,Nt); U(1:(5/dt))=0.5; 
        else
            U = 0.075*ones(1,Nt); U((20/dt):(30/dt))=0;
        end

        %State equations
        State = @(t,y,U) [(ba+pa*y(1)*(1-y(1)-y(2)-y(3))-ma*y(1)-mau*U*y(1)); 
            (pp*y(2)*(1-y(1)-y(2)-y(3)) -dp*y(2) +dn*y(3) -dpu*U*y(2) -mp*y(2)-mpu*U*y(2) - alpha_*y(2)/(gamma_ + y(2)+y(3)) ); 
            (pn*y(3)*(1-y(1)-y(2)-y(3)) +dp*y(2) -dn*y(3) +dpu*U*y(2) -mn*y(3)) - alpha_*y(3)/(gamma_ + y(2)+y(3)) ];
        
        %Forward sweep using fourth-order Runge-Kutta scheme
        i = 0; %Initialise loop variable
        t = 0; %Initialise time for forward sweep
        while i < Nt-1
            t = t + dt;
            i=i+1;
            k1 = State(t,y(:,i),U(i));
            k2 = State(t,y(:,i)+dt*k1/2,0.5*(U(i)+U(i+1)));
            k3 = State(t,y(:,i)+dt*k2/2,0.5*(U(i)+U(i+1)));
            k4 = State(t,y(:,i)+dt*k3,U(i+1));
            y(:,i+1) = y(:,i) + (dt/6)*(k1+2*k2+2*k3+k4);
        end

        figure
        set(gca, 'ColorOrder', colours);
        hold on
        box on
        line1 = plot(t_y,y(1,:),'LineWidth',2);
        line2 = plot(t_y,y(2,:),'LineWidth',2);
        line3 = plot(t_y,y(3,:),'LineWidth',2);
        line4 = plot(t_y,U,'k','LineWidth',2);
        hL = legend([line1,line2,line3,line4],{'A','P','N','u'},'Location','northeast');
        ylabel('State','fontsize',18);
        xlabel('Time','fontsize',18);
        axis([0,Tfinal,0,1])
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 18)
        set(gcf, 'Position',  [100, 100, 400, 480]);

        saveas(gcf,append(nm,'.png') );
        close(gcf);

    end
end






