function dPdt = model(t,P,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope)
% Zhuojun Yu: This function is the ODE system of CPG-FB model

% CPG parameters
Iapp=0.8;                      
Vk=-80; Vl=-50; Vca=100;       
gk=0.02; gl=0.005; gca=0.015; 
c=1;                           
E1=0; E3=0;    
E2=15; E4=15;
Eslope=2;                     
phi=0.0005;          
Esyn=-80;        
% Muscle parameters             
tau=2.45;                   
beta=0.703;             
a0=0.165;               
g=2;                    
F0=10;              
b=4000;    

V1=P(1); V2=P(2); N1=P(3); N2=P(4); A1=P(5); A2=P(6); x=P(7); y=P(8);

minf1 = 0.5*(1+tanh((V1-E1)/E2));
minf2 = 0.5*(1+tanh((V2-E1)/E2));
winf1 = 0.5*(1+tanh((V1-E3)/E4));
winf2 = 0.5*(1+tanh((V2-E3)/E4));
tauw1 = 1/cosh((V1-E3)/(2*E4));
tauw2 = 1/cosh((V2-E3)/(2*E4));

sinffw1 = 0.5*(1+tanh((V1-Ethresh)/Eslope));
sinffw2 = 0.5*(1+tanh((V2-Ethresh)/Eslope));
Isyn1 = gsyn*sinffw2*(V1-Esyn);    
Isyn2 = gsyn*sinffw1*(V2-Esyn);

L1 = 10+x; L2 = 10-x;

u1 = (1/2)*V1; u2 = (1/2)*V2;
U1 = (1.03-4.31*exp(-0.198*u1))*(u1>=8);
U2 = (1.03-4.31*exp(-0.198*u2))*(u2>=8);
LT1 = -3*sqrt(3)*(L1-1)*(L1-5)*(L1-15)/2/625;
LT2 = -3*sqrt(3)*(L2-1)*(L2-5)*(L2-15)/2/625;

a1 = g*(A1-a0); a2 = g*(A2-a0);
F1 = (F0*a1*LT1)*(u1>=8)*(a1>=0);
F2 = (F0*a2*LT2)*(u2>=8)*(a2>=0);

sinffb1 = 0.5*(1-tanh((L1-L0)/Lslope)); 
sinffb2 = 0.5*(1-tanh((L2-L0)/Lslope)); 
Ifb1 = gfb*sinffb2*(V1-Efb);   
Ifb2 = gfb*sinffb1*(V2-Efb);  

dV1dt = (Iapp-gca*minf1*(V1-Vca)-gk*N1*(V1-Vk)-gl*(V1-Vl)-Isyn1-Ifb1)/c;
dV2dt = (Iapp-gca*minf2*(V2-Vca)-gk*N2*(V2-Vk)-gl*(V2-Vl)-Isyn2-Ifb2)/c;
dN1dt = phi*(winf1-N1)/tauw1;
dN2dt = phi*(winf2-N2)/tauw2;
dA1dt = (1/tau)*(U1-(beta+(1-beta)*U1)*A1);
dA2dt = (1/tau)*(U2-(beta+(1-beta)*U2)*A2);
if V1 > Ethresh
    r=1;
    dxdt = (1/b)*(F2-F1+kappa*F_ell);
else
    r=0;
    dxdt = (1/b)*(F2-F1);
end
dydt = -r*dxdt;

dPdt = [dV1dt;dV2dt;dN1dt;dN2dt;dA1dt;dA2dt;dxdt;dydt];
end