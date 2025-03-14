function [cost,U,A,P,N] = steady_state_cost(params,CostFn,controlLevel, tol)

%steady_state_cost: For the given control level, this 
% function calculates the infinite time horizon time-averaged cost.
% If there are more than 1 stable cancer states associated with the given 
% control level, assumes cancer has been returned to the lowest.
% returns this together with the associated U,A,P,N
% First check if the N=P=0 solution is stable for the given control level
% If so, return based on that.
% Otherwise, if alpha=0 then calculate the exact steady state P,N>0
% solution (must be stable) and return based on that.
% Otherwise, numerically calculate the P,N>0 steady state using 
% steadystates_dara_immune function, and return based on that.

% controlLevel = 0.1 %test

U=controlLevel;
for key = params.keys()
    eval(append(key{1}," = params('",key{1},"');"));
end
% adjust the parameters to incorporate the control and the immune response
% at P,N ~ 0
ma = ma + mau*controlLevel;
mp = mp + mpu*controlLevel + alpha_/gamma_;
dp = dp + dpu*controlLevel;
mn = mn +alpha_/gamma_;

% C=1-A-P-N
% C_bound is maximum equilibrium C, corresponding to N=P=0
% Cc is equilibrium C value for P,N>0: physically realisable if Cc<C_bound.
% With parameters adjusted to include immune response at N,P ~0, Cc >= C_bound
% becomes a test for stability of the N=P=0 equilibrium (with immune
% response and control)

C_bound = (pa+ma)/(2*pa) - sqrt(((pa-ma)/(2*pa))^2 + ba/pa);

%if (dp==0) && (dn==0) % P and N can't coexist; assume most fit survives (for stable equilib)

if pn==0
    Cc = (dp+mp-dp*dn/(dn+mn))/pp;
elseif pp==0
    Cc = (dn+mn-dp*dn/(dp+mp))/pn;
else
    Cc = (dp+mp)/(2*pp) + (dn+mn)/(2*pn) - sqrt( ((dp+mp)/(2*pp) - (dn+mn)/(2*pn))^2 + dp*dn/(pp*pn));
    % if (dp==0) && (dn==0), gives Cc = min(mp/pp,mn/pn) [correct]
    % if (dp>0) && (dn==0), gives Cc = min((dp+mp)/pp,mn/pn) [correct]
end

if Cc >= C_bound
    A = 1-C_bound;
    P=0;
    N=0;
elseif alpha_==0
    %"ctrl: "+
    % controlLevel
    A = ba/(ma-pa*Cc);
    PNsum=1-Cc-A;
    fN=dp;
    fP=dn+mn-pn*Cc;
    % do dn,dp=0 cases explicitly to avoid any traps
    if (dp==0)
        if (dn==0)
            if (mp/pp > mn/pn)
                fN=1;
                fP=0;
            else
                fN=0;
                fP=1;
            end
        else

        end
    elseif (dn==0)
        if ((dp+mp)/pp > mn/pn)
            fP=0;
            fN=1;
        else
            % in this case original values work
            fN=dp;
            fP=mn-pn*Cc;
        end

    end

    % [dn , dp+mp-pp*Cc]
    P=PNsum*fP/(fP+fN);
    N=PNsum*fN/(fP+fN);
else
    APN_vals = steadystates_dara_immune(params,controlLevel,tol); 
    A = APN_vals{end,'A'};
    P = APN_vals{end,'P'};
    N = APN_vals{end,'N'};
end
cost=CostFn(controlLevel,P+N);

end
