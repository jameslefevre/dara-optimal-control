function [APN_vals,C_bound,Cc] = steadystates_dara_immune(params,controlLevel,tol)
%steadystates_dara_immune Calculates the steady state(s) of the dara/MM
%system with immune response
%   Steady state is calculated given a constant control value u
%   Set to 0 to get SS sans control
% Initially get the steady state without immune response.
% (If the coexist SS exists with immune response, then it exists without)
% Any coexist SS is then run through simulation until change is <tol in all
% 3 state variables

% export parameters from params map to local namespace
for key = params.keys()
    eval(append(key{1}," = params('",key{1},"');"));
end
if controlLevel>0
    ma = ma + mau*controlLevel;
    mp = mp + mpu*controlLevel; % have used mpu=1, and don't always mention
    dp = dp + dpu*controlLevel;
end
% leave t out of the state equation, since SS assumes lack of time dependence 
State = @(y,U) [(ba+pa*y(1)*(1-y(1)-y(2)-y(3))-ma*y(1)-mau*U*y(1)); 
    (pp*y(2)*(1-y(1)-y(2)-y(3)) -dp*y(2) +dn*y(3) -dpu*U*y(2) -mp*y(2)-mpu*U*y(2) - alpha_*y(2)/(gamma_ + y(2)+y(3)) ); 
    (pn*y(3)*(1-y(1)-y(2)-y(3)) +dp*y(2) -dn*y(3) +dpu*U*y(2) -mn*y(3)) - alpha_*y(3)/(gamma_ + y(2)+y(3)) ];


C_bound = (pa+ma)/(2*pa) - sqrt(((pa-ma)/(2*pa))^2 + ba/pa);
A_nc = 1-C_bound;
if pn==0
    Cc = (dp+mp-dp*dn/(dn+mn))/pp;
elseif pp==0
    Cc = (dn+mn-dp*dn/(dp+mp))/pn;
else
    Cc = (dp+mp)/(2*pp) + (dn+mn)/(2*pn) - sqrt( ((dp+mp)/(2*pp) - (dn+mn)/(2*pn))^2 + dp*dn/(pp*pn));
end

sol_num = 1;
if Cc<C_bound
    sol_num=2;
end
APN_vals = table('Size',[sol_num,3],'VariableTypes',["double","double","double"] , 'VariableNames',["A","P","N"] );
APN_vals{1,"A"}=A_nc;
APN_vals{1,"P"}=0;
APN_vals{1,"N"}=0;
if Cc<C_bound
    
    A = ba/(ma-pa*Cc);
    PNsum=1-Cc-A;
    fN=dp;
    fP=dn+mn-pn*Cc;
    y = [A,PNsum*fP/(fP+fN),PNsum*fN/(fP+fN)];
    for dummy = 1:1e6
        k1 = State(y,0);
        k2 = State(y+dt*k1/2,0);
        k3 = State(y+dt*k2/2,0);
        k4 = State(y+dt*k3,0);
        diff = (dt/6)*(k1+2*k2+2*k3+k4);
        y = y + (dt/6)*(k1+2*k2+2*k3+k4);
        if max(abs(diff))<tol
            break;
        end
    end
    APN_vals{2,"A"}=y(1);
    APN_vals{2,"P"}=y(2);
    APN_vals{2,"N"}=y(3);
end
end