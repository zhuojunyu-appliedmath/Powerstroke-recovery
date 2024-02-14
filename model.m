function dPdt = model(t,P,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf)

% Neuron numbering: 1, RG-F; 2, RG-E; 3, In-F; 4, In-E; 5, PF-F; 6, PF-E; 
% 7, Int; 8, Inab-E; 9, Mn-F; 10, Mn-E; 
% Feedback numbering: 1, Iaf; 2, IIf; 3, Iae; 4, Ibe

V1=P(1); V2=P(2); V3=P(3); V4=P(4); V5=P(5); V6=P(6); V7=P(7); V8=P(8); V9=P(9); V10=P(10);
h1=P(11); h2=P(12); h5=P(13); h6=P(14); h9=P(15); h10=P(16); 
q=P(17); v=P(18); 

q=max(q,0); q=min(q,pi); %angle restrictions

%% Parameters
%CPG
C=20;
Ena=55; Ek=-80; Esyne=-10; Esyni=-70; 
El1=-64;    %for RG, PF and Mn
El2=-60;    %otherwise
gk=4.5; gl=1.6; gsyne=10; gsyni=10; 
gnap1=3.5;  %for RG
gnap2=0.5;  %for PF
gnap3=0.3;  %for Mn
Vhalf=-30; Vth=-50;
k1=3;       %for Mn
k2=8;       %otherwise
d1=1.4;       %constant supra-spinal drive
c1=0.08; c2=0.08; c5=0.4; c6=0.4; %supra-spinal drive strengths
d2=0.18;    %external drive to int
%excitatory connection 
a13=0.41; a24=0.41;
a15=0.7; a26=0.7;
a59=1.95;
a610=1.3; a810=0.82;
a68=0.35;
%inhibitory connection 
b41=2.2; b32=2.2;
b45=6.6; b36=6.6;
b47=2.8; 
b78=0.55;
%feedback connection
w11=0.06; w21=0.0348; 
w32=0.06; w42=0.066;
w13=0.27; w23=0.1566;
w34=0.44; w44=0.484;
w15=0.19; w25=0.1102;
w36=0.1; w46=0.11;
w38=0.16; w48=0.176;
%limb
m=300; ls=300; g=0.00981;
I=m*ls^2/3;
K=0.5*m*g*ls;
b=0.002;
MGRmax=585;
%muscles
a1=60; a2=7;
Lopt=68;
Ffmax=-72.5; Femax=-37.7;
%feedback
Lth_Ia=60.007; 
Lth_II=58.457;  %minimal muscle length evoking afferent activation
Fth=3.393;
rhov=0.6;
kvIa=6.2; kdIa=2;
kEMGIa=0.06; kEMGII=0.06;
CIa=0.026; CII=0;
kIb=1; kdII=1.5;

%% Functions
%(in)activation
mnap = @(V) (1+exp(-(V+47.1)/3.1))^(-1);
mk = @(V) (1+exp(-(V+44.5)/5))^(-1);
hinf = @(V) (1+exp((V+51)/4))^(-1);
tauh = @(V) 600*(cosh((V+51)/8))^(-1);
%currents
Inap = @(V,h,gnap) gnap*mnap(V)*h*(V-Ena);
Ik = @(V) gk*mk(V)^4*(V-Ek);
Il = @(V,El) gl*(V-El);
f = zeros(10,1);  %neuron outputs
for i=1:8
    f(i) = (P(i)>=Vth)*(1+exp(-(P(i)-Vhalf)/k2))^(-1);
end
for i=9:10
    f(i) = (P(i)>=Vth)*(1+exp(-(P(i)-Vhalf)/k1))^(-1);
end
%muscle 
Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q));    %flexor muscle length
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); %extensor muscle length
lf = Lf/Lopt; le = Le/Lopt; %normalized muscle lengths
hf = a1*a2*sin(q)/Lf;       %flexor muscle moment arm
he = a1*a2*sin(pi-q)/Le;    %extensor muscle moment arm
vmf = v*hf; vme = -v*he;    %muslce velocities
Fl = @(l) exp(-(abs((l^2.3-1)/1.26))^1.62); %force dependence on muscle length
Fv = @(vm,l) (vm<0)*(-0.69-0.17*vm)/(vm-0.69)+(vm>=0)*(0.18-(-5.34*l^2+8.41*l-4.7)*vm)/(vm+0.18); %force dependence on velocity
Fp = @(l) 3.5*log(exp((l-1.4)/0.05)+1)-0.02*(exp(-18.7*(l-0.79))-1);  %passive force
Ff = Ffmax*(f(9)*Fl(lf)*Fv(vmf,lf)+Fp(lf));  %flexor muscle force
Fe = Femax*(f(10)*Fl(le)*Fv(vme,le)+Fp(le)); %extensor muscle force
Mf = Ff*hf; Me = -Fe*he;    %muslce moments
MGR = -MGRmax*cos(q-kappa)*(v>0);  %moment of ground reaction force
%feedback currents
Iaf = s_Iaf*max(sign(vmf)*kvIa*abs(vmf/Lth_Ia)^rhov+kdIa*max((Lf-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(9)+CIa,0);
Iae = s_Iae*max(sign(vme)*kvIa*abs(vme/Lth_Ia)^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(10)+CIa,0);
%Iaf = s_Iaf*max(sign(vmf)*kvIa*abs(vmf/Lth_Ia)^rhov+kdIa*(1+tanh((Lf-Lth_Ia)/Lslope))/20+kEMGIa*f(9)+CIa,0);
%Iae = s_Iae*max(sign(vme)*kvIa*abs(vme/Lth_Ia)^rhov+kdIa*(1+tanh((Le-Lth_Ia)/Lslope))/20+kEMGIa*f(10)+CIa,0);
%Iaf = s_Iaf*kdIa*(1+tanh((Lf-Lth_Ia)/Lslope))/20;
%Iae = s_Iae*kdIa*(1+tanh((Le-Lth_Ia)/Lslope))/20;
Ibe = -s_Ibe*kIb/Femax*max(-Fe-Fth,0);%kIb*(Fe<=-Fth)*(-Fe-Fth)/(-Femax);
IIf = s_IIf*(kdII*max((Lf-Lth_II)/Lth_II,0)+kEMGII*f(9)+CII);
%excitatory and inhibitory inputs
Isyne1 = gsyne*(V1-Esyne)*(c1*d1+w11*Iaf+w21*IIf);
Isyni1 = gsyni*(V1-Esyni)*b41*f(4);
Isyne2 = gsyne*(V2-Esyne)*(c2*d1+w32*Iae+w42*Ibe);
Isyni2 = gsyni*(V2-Esyni)*b32*f(3);
Isyne3 = gsyne*(V3-Esyne)*(a13*f(1)+w13*Iaf+w23*IIf);
Isyne4 = gsyne*(V4-Esyne)*(a24*f(2)+w34*Iae+w44*Ibe);
Isyne5 = gsyne*(V5-Esyne)*(c5*d1+a15*f(1)+w15*Iaf+w25*IIf);
Isyni5 = gsyni*(V5-Esyni)*b45*f(4);
Isyne6 = gsyne*(V6-Esyne)*(c6*d1+a26*f(2)+w36*Iae+w46*Ibe);
Isyni6 = gsyni*(V6-Esyni)*b36*f(3);
Isyne7 = gsyne*(V7-Esyne)*d2;
Isyni7 = gsyni*(V7-Esyni)*b47*f(4);
Isyne8 = gsyne*(V8-Esyne)*(a68*f(6)+w38*Iae+w48*Ibe);
Isyni8 = gsyni*(V8-Esyni)*b78*f(7);
Isyne9 = gsyne*(V9-Esyne)*a59*f(5);
Isyne10 = gsyne*(V10-Esyne)*(a610*f(6)+a810*f(8));    
    
%% ODE system
dV1dt = (-Inap(V1,h1,gnap1)-Ik(V1)-Il(V1,El1)-Isyne1-Isyni1)/C;
dV2dt = (-Inap(V2,h2,gnap1)-Ik(V2)-Il(V2,El1)-Isyne2-Isyni2)/C;
dV3dt = (-Il(V3,El2)-Isyne3)/C;
dV4dt = (-Il(V4,El2)-Isyne4)/C;
dV5dt = (-Inap(V5,h5,gnap2)-Ik(V5)-Il(V5,El1)-Isyne5-Isyni5)/C;
dV6dt = (-Inap(V6,h6,gnap2)-Ik(V6)-Il(V6,El1)-Isyne6-Isyni6)/C;
dV7dt = (-Il(V7,El2)-Isyne7-Isyni7)/C;
dV8dt = (-Il(V8,El2)-Isyne8-Isyni8)/C;
dV9dt = (-Inap(V9,h9,gnap3)-Ik(V9)-Il(V9,El1)-Isyne9)/C;
dV10dt = (-Inap(V10,h10,gnap3)-Ik(V10)-Il(V10,El1)-Isyne10)/C;
dh1dt = (hinf(V1)-h1)/tauh(V1);
dh2dt = (hinf(V2)-h2)/tauh(V2);
dh5dt = (hinf(V5)-h5)/tauh(V5);
dh6dt = (hinf(V6)-h6)/tauh(V6);
dh9dt = (hinf(V9)-h9)/tauh(V9);
dh10dt = (hinf(V10)-h10)/tauh(V10);
dqdt = v;
dvdt = (K*cos(q)+Mf+Me+MGR)/I-b*v;

dPdt = [dV1dt;dV2dt;dV3dt;dV4dt;dV5dt;dV6dt;dV7dt;dV8dt;dV9dt;dV10dt;...
        dh1dt;dh2dt;dh5dt;dh6dt;dh9dt;dh10dt;...
        dqdt;dvdt];

end
