% find period, power-stroke duration, and initial and end point on power stroke and recovery for a limit cycle
function [T0,T0_ps,init_ps,end_ps,init_re,end_re] = phases(gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,init)

tF1 = 15000; tF2 = 2800;

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
options_ps = odeset('Events',@event_ps,'RelTol',1e-8,'AbsTol',1e-8);
options_re = odeset('Events',@event_re,'RelTol',1e-8,'AbsTol',1e-8);

[~,P] = ode15s(@model,[0 tF1],init,[],gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
[~,~,~,init_ps,~] = ode15s(@model,[0 tF2],P(end,:),options_ps,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
[~,~,T0,init_ps,~] = ode15s(@model,[0 tF2],init_ps(end,:),options_ps,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
[~,~,T0_ps,init_re,~] = ode15s(@model,[0 tF2],init_ps(1,:),options_re,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);

init_ps = init_ps(1,1:7);
init_re = init_re(1,1:7);
T0_ps=T0_ps(1);

dt=0.01;
[~,P] = ode15s(@model,[0:dt:T0_ps],[init_ps 0],options,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
end_ps = P(end,1:7);
[~,P] = ode15s(@model,[0:dt:T0],[init_ps 0],options,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
end_re = P(end,1:7);

end
