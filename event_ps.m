function [value,isterminal,direction] = event_ps(t,P,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope)
% This function finds the initial point on power-stroke phase

V1=P(1); 

value=V1-Ethresh;
isterminal=0;
direction=1;

end