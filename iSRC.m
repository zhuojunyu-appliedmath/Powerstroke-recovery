function dPdt = iSRC(t,P,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,nu1_ps,nu1_re)
% This function is used to compute (piecewise) iSRC

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
% Muscle model parameters             
tau=2.45;                   
beta=0.703;             
a0=0.165;               
g=2;                    
F0=10;              
b=4000;         

V1=P(1); V2=P(2); N1=P(3); N2=P(4); A1=P(5); A2=P(6); x=P(7); y=P(8);
src1=P(9); src2=P(10); src3=P(11); src4=P(12); src5=P(13); src6=P(14); src7=P(15);
src=[src1;src2;src3;src4;src5;src6;src7];

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

sinffb1 = 0.5*(1-tanh((L1-L0)/Lslope)); 
sinffb2 = 0.5*(1-tanh((L2-L0)/Lslope)); 
Ifb1 = gfb*sinffb2*(V1-Efb);   
Ifb2 = gfb*sinffb1*(V2-Efb);  

u1 = (1/2)*V1; u2 = (1/2)*V2;
U1 = (1.03-4.31*exp(-0.198*u1))*(u1>=8);
U2 = (1.03-4.31*exp(-0.198*u2))*(u2>=8);
LT1 = -3*sqrt(3)*(L1-1)*(L1-5)*(L1-15)/2/625;
LT2 = -3*sqrt(3)*(L2-1)*(L2-5)*(L2-15)/2/625;

a1 = g*(A1-a0); 
F1 = (F0*a1*LT1)*(u1>=8)*(a1>=0);
a2 = g*(A2-a0);
F2 = (F0*a2*LT2)*(u2>=8)*(a2>=0);

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

% iSRC

dminf1 = sech((V1-E1)/E2)^2/2/E2;
dminf2 = sech((V2-E1)/E2)^2/2/E2;
dwinf1 = sech((V1-E3)/E4)^2/2/E4;
dwinf2 = sech((V2-E3)/E4)^2/2/E4;

dsinffw1 = sech((V1-Ethresh)/Eslope)^2/2/Eslope;
dsinffw2 = sech((V2-Ethresh)/Eslope)^2/2/Eslope;

dsinffb1 = -sech((L1-L0)/Lslope)^2/2/Lslope;
dsinffb2 = -sech((L2-L0)/Lslope)^2/2/Lslope;

dlambda1 = phi*sinh((V1-E3)/E4)/2/E4;
dlambda2 = phi*sinh((V2-E3)/E4)/2/E4;

dLT1 = -3*sqrt(3)*(3*L1^2-42*L1+95)/2/625;
dLT2 = -3*sqrt(3)*(3*L2^2-42*L2+95)/2/625;

if u1>=8
    dU1 = 4.31*exp(-0.198*u1)*0.099;
    if a1>=0
        dF1dA1 = F0*g*LT1;
        dF1dx = F0*a1*dLT1*0.8;
    else 
        dF1dA1 = 0;
        dF1dx = 0;
    end
else
    dU1 = 0;
    dF1dA1 = 0;
    dF1dx = 0;
end

if u2>=8
    dU2 = 4.31*exp(-0.198*u2)*0.099;
    if a2>=0
        dF2dA2 = F0*g*LT2;
        dF2dx = F0*a2*dLT2*(-0.8);
    else 
        dF2dA2 = 0;
        dF2dx = 0;
    end
else
    dU2 = 0;
    dF2dA2 = 0;
    dF2dx = 0;
end

% Jacobian
Jac = [(-gl-gca*dminf1*(V1-Vca)-gca*minf1-gk*N1-gsyn*sinffw2-gfb*sinffb2)/c  -gsyn*dsinffw2*(V1-Esyn)/c  -gk*(V1-Vk)/c  0  0  0  gfb*dsinffb2*(V1-Efb)/c;...
       -gsyn*dsinffw1*(V2-Esyn)/c  (-gl-gca*dminf2*(V2-Vca)-gca*minf2-gk*N2-gsyn*sinffw1-gfb*sinffb1)/c  0  -gk*(V2-Vk)/c  0  0  -gfb*dsinffb1*(V2-Efb)/c;...
       dlambda1*(winf1-N1)+phi*dwinf1/tauw1  0  -phi/tauw1  0  0  0  0;...
       0  dlambda2*(winf2-N2)+phi*dwinf2/tauw2  0  -phi/tauw2  0  0  0;...
       (dU1-(1-beta)*dU1*A1)/tau  0  0  0  -(beta+(1-beta)*U1)/tau  0  0;...
       0  (dU2-(1-beta)*dU2*A2)/tau  0  0  0  -(beta+(1-beta)*U2)/tau  0;...
       0  0  0  0  -dF1dA1/b  dF2dA2/b  (dF2dx-dF1dx)/b];

% Vector field    
VF = [dV1dt;dV2dt;dN1dt;dN2dt;dA1dt;dA2dt;dxdt];

if r==1
    dVF = [0;0;0;0;0;0;F_ell/b];
    dsrcdt = Jac*src+nu1_ps*VF+dVF;
else
    dsrcdt = Jac*src+nu1_re*VF;
end

dPdt = [dV1dt;dV2dt;dN1dt;dN2dt;dA1dt;dA2dt;dxdt;dydt;dsrcdt];

end