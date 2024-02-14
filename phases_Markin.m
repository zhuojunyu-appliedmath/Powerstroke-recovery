% find period, stance duration, and initial and end point on stance and swing for a limit cycle
function [T0,T0_st,init_st,end_st,init_sw,end_sw] = phases(init,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf)

tF1 = 9000; tF2 = 1800;

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
options_st = odeset('Events',@sw_to_st,'RelTol',1e-8,'AbsTol',1e-8);
options_sw = odeset('Events',@st_to_sw,'RelTol',1e-8,'AbsTol',1e-8);
[~,P] = ode15s(@model,[0 tF1],init,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
[~,~,~,init_st,~] = ode15s(@model,[0 tF2],P(end,:),options_st,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
[~,~,T0,init_st,~] = ode15s(@model,[0 tF2],init_st(end,:),options_st,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
[~,~,T0_st,init_sw,~] = ode15s(@model,[0 tF2],init_st(1,:),options_sw,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);

init_st = init_st(1,:);
init_sw = init_sw(1,:);
T0_st=T0_st(1);

dt=0.1;
[~,P] = ode15s(@model,[0:dt:T0_st],init_st,options,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
end_st = P(end,:);
[~,P] = ode15s(@model,[0:dt:T0],init_st,options,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
end_sw = P(end,:);

end
