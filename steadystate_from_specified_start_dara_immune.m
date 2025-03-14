function [A,P,N] = steadystate_from_specified_start_dara_immune(params,A0,P0,N0,controlLevel,tol)
% steadystate_from_specified_start_dara_immune runs a simulation for a
% given start state and specified fixed control level, until stable within
% specified tolerance.
% Returns final state

% export parameters from params map to local namespace
for key = params.keys()
    eval(append(key{1}," = params('",key{1},"');"));
end
if controlLevel>0
    ma = ma + mau*controlLevel;
    mp = mp + mpu*controlLevel; % have used mpu=1, and don't always mention
    dp = dp + dpu*controlLevel;
end

State = @(y,U) [(ba+pa*y(1)*(1-y(1)-y(2)-y(3))-ma*y(1)-mau*U*y(1)); 
    (pp*y(2)*(1-y(1)-y(2)-y(3)) -dp*y(2) +dn*y(3) -dpu*U*y(2) -mp*y(2)-mpu*U*y(2) - alpha_*y(2)/(gamma_ + y(2)+y(3)) ); 
    (pn*y(3)*(1-y(1)-y(2)-y(3)) +dp*y(2) -dn*y(3) +dpu*U*y(2) -mn*y(3)) - alpha_*y(3)/(gamma_ + y(2)+y(3)) ];


y = [A0;P0;N0];
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
A=y(1);
P=y(2);
N=y(3);