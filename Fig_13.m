% Fig. 13: nullcline configuration in the default system

clear; clc; close all

%% Parameters

% Parameters
gnap1=3.5; gk=4.5; gl=1.6; gsyne=10; gsyni=10; 
Ena=55; El1=-64; Ek=-80; Esyne=-10; Esyni=-70; 
c1=0.08; c2=0.08;
d1=1.4;      
Vhalf=-30; Vth=-50;
k1=3; k2=8;      
w11=0.06; w21=0.0348; 
w32=0.06; w42=0.066;
w13=0.27; w23=0.1566;
w34=0.44; w44=0.484;
b41=2.2; b32=2.2;
a1=60; a2=7;
Lopt=68;
Femax=-37.7;
Lth_Ia=60.007; Lth_II=58.457;  
Fth=3.393;
rhov=0.6;
kvIa=6.2; kdIa=2;
kEMGIa=0.06; kEMGII=0.06;
CIa=0.026; CII=0;
kIb=1; kdII=1.5;

maxV=-30; minV=-65;  
V = minV:0.1:maxV;

s_Iaf=1;
s_Iae=1;
s_Ibe=1;
s_IIf=1;

kappa=0;

% Functions
mnap = @(V) (1+exp(-(V+47.1)/3.1)).^(-1);
mk = @(V) (1+exp(-(V+44.5)/5)).^(-1);
f3 = @(V) (V>=Vth).*(1+exp(-(V-Vhalf)/k2)).^(-1);
f4 = @(V) (V>=Vth).*(1+exp(-(V-Vhalf)/k2)).^(-1);
f9 = @(V) (V>=Vth).*(1+exp(-(V-Vhalf)/k1)).^(-1);
f10 = @(V) (V>=Vth).*(1+exp(-(V-Vhalf)/k1)).^(-1);
Ik = @(V) gk*mk(V).^4.*(V-Ek);
Il = @(V,El) gl*(V-El);
Fl = @(l) exp(-(abs((l.^2.3-1)/1.26)).^1.62); 
Fv = @(vm,l) (vm<0).*(-0.69-0.17*vm)./(vm-0.69)+(vm>=0).*(0.18-(-5.34*l.^2+8.41*l-4.7).*vm)./(vm+0.18); 
Fp = @(l) 3.5*log(exp((l-1.4)/0.05)+1)-0.02*(exp(-18.7*(l-0.79))-1); 
h_null = (1+exp((V+51)/4)).^(-1); %h-nullcline


%% Solution
init_st=[-64.8809378117064,-36.9594157966064,-58.6472908751815,-27.7472995924213,-62.2355322784133,-33.0647767211889,-63.6072521709523,-31.9808006359821,-63.9624843683281,-32.3185145652701,0.533520879043655,0.301667961622656,0.392933844976166,0.167281647228994,0.403832124325481,0.232029016617407,1.29920012370126,8.29921889103028e-19];
T0=1035.22; T0_st=719.03;
[tspan,P] = ode15s(@model,[0:.1:T0],init_st,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
V1=P(:,1); V2=P(:,2); V3=P(:,3); V4=P(:,4); V9=P(:,9); V10=P(:,10); 
h1=P(:,11); h2=P(:,12); 
q=P(:,17); v=P(:,18);


%% Nullclines for RG

%special points as In-E jumps up across Vth
for k=1:length(V4)-1
    if V4(k)<V4(k+1) && V4(k)-Vth<-0.7 && V4(k)-Vth>-0.8
        ind1=k; 
    end
end
for k=1:length(V4)-1
    if V4(k)<V4(k+1) && V4(k)-Vth>0.7 && V4(k)-Vth<0.8
        ind2=k;
    end
end

%feedback
Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); 
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
le = Le/Lopt; 
vmf = v.*hf; vme = -v.*he;    
Fe = Femax*(f10(V10).*Fl(le).*Fv(vme,le)+Fp(le)); 

Iaf = s_Iaf*max(sign(vmf)*kvIa.*abs(vmf/Lth_Ia).^rhov+kdIa*max((Lf-Lth_Ia)/Lth_Ia,0)+kEMGIa*f9(V9)+CIa,0);
IIf = s_IIf*(kdII*max((Lf-Lth_II)/Lth_II,0)+kEMGII*f9(V9)+CII);
Iae = s_Iae*max(sign(vme)*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f10(V10)+CIa,0);
Ibe = -s_Ibe*kIb/Femax*max(-Fe-Fth,0);

% RG-F
Isyne1_a1 = gsyne*(V-Esyne)*(c1*d1+w11*Iaf(ind1)+w21*IIf(ind1));  %at active state (before In-F jumping)
Isyne1_a2 = gsyne*(V-Esyne)*(c1*d1+w11*Iaf(ind2)+w21*IIf(ind2));  %at active state (after In-F jumping)
Isyni2_a2 = gsyni*(V-Esyni)*b41*f4(V4(ind2));

V1_a_null_l = (-Ik(V)-Il(V,El1)-Isyne1_a1)./(gnap1*mnap(V).*(V-Ena)); %active V-nullcline (before In-F jumping)
V1_a_null_2 = (-Ik(V)-Il(V,El1)-Isyne1_a2-Isyni2_a2)./(gnap1*mnap(V).*(V-Ena)); %active V-nullcline (after In-F jumping)

% RG-E
Isyne2_s1 = gsyne*(V-Esyne)*(c2*d1+w32*Iae(ind1)+w42*Ibe(ind1)); %at silent state (before In-F jumping)
Isyne2_s2 = gsyne*(V-Esyne)*(c2*d1+w32*Iae(ind2)+w42*Ibe(ind2)); %at silent state (after In-F jumping)

Isyni2_s1 = gsyni*(V-Esyni)*b32*f3(V3(ind1)); %fully inhibitd by In-F
Isyni2_s2 = gsyni*(V-Esyni)*b32*f3(V3(ind2)); 

V2_s_null_1 = (-Ik(V)-Il(V,El1)-Isyne2_s1-Isyni2_s1)./(gnap1*mnap(V).*(V-Ena)); %silent V-nullclines (before In-F jumping)
V2_s_null_2 = (-Ik(V)-Il(V,El1)-Isyne2_s2-Isyni2_s2)./(gnap1*mnap(V).*(V-Ena)); %silent V-nullclines (after In-F jumping)


%% Plot nullclines

figure(1)
 
subplot(2,2,1)
plot(V,h_null,'-g','LineWidth',1.5); hold on
plot(V,V1_a_null_l,'-b','LineWidth',1.5);
plot(V,V2_s_null_1,'-r','LineWidth',1.5);
plot(V1(ind1),h1(ind1),'.b','MarkerSize',25); 
plot(V2(ind1),h2(ind1),'.r','MarkerSize',25); 
plot(Vth*ones(length([0:.1:1]),1),[0:.1:1],'-m','LineWidth',1.5); 
axis([minV,maxV,0,1])
xlabel('V_{RG}'); ylabel('h_{RG}')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(V3(ind1),0,'*b','MarkerSize',7); hold on
plot(V4(ind1),0,'*r','MarkerSize',7); 
plot(Vth*ones(length([-0.1:.01:.1]),1),[-0.1:.01:.1],'-m','LineWidth',1.5); 
axis([minV,maxV,-.1,.1])
xlabel('V_{In}'); 
set(gca,'FontSize',13)

subplot(2,2,2)
plot(V,h_null,'-g','LineWidth',1.5); hold on
plot(V,V1_a_null_2,'-b','LineWidth',1.5);
plot(V,V2_s_null_2,'-r','LineWidth',1.5);
plot(V1(ind2),h1(ind2),'.b','MarkerSize',25); 
plot(V2(ind2),h2(ind2),'.r','MarkerSize',25); 
plot(Vth*ones(length([0:.1:1]),1),[0:.1:1],'-m','LineWidth',1.5); 
axis([minV,maxV,0,1])
xlabel('V_{RG}'); ylabel('h_{RG}')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(V3(ind2),0,'ob','MarkerSize',7); hold on
plot(V4(ind2),0,'or','MarkerSize',7); 
plot(Vth*ones(length([-0.1:.01:.1]),1),[-0.1:.01:.1],'-m','LineWidth',1.5); 
axis([minV,maxV,-.1,.1])
xlabel('V_{In}'); 
set(gca,'FontSize',13)
