% Fig. 10

% Neuron numbering: 1, RG-F; 2, RG-E; 3, In-F; 4, In-E; 5, PF-F; 6, PF-E; 
% 7, Int; 8, Inab-E; 9, Mn-F; 10, Mn-E; 

clear; clc; close all;

%feedback strength
s_Iaf=1;
s_Iae=1;
s_Ibe=1;
s_IIf=1;

kappa=0;
init_st=[-64.8809361683031,-36.9594183760335,-58.6472831890731,-27.7473050183001,-62.2355295469095,-33.0647778390756,-63.6072508638440,-31.9808035411915,-63.9624843671704,-32.3185157661491,0.533521095965132,0.301667837348812,0.392933826237958,0.167281516515830,0.403832356834355,0.232028801749969,1.29920016313509,4.57591682652905e-17];
T0=1035.2913; T0_st=719.0293; tF=2*T0;

%% solve ODEs
[tspan,P] = ode15s(@model,[0 tF],init_st,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
V1=P(:,1); V2=P(:,2); V3=P(:,3); V4=P(:,4); V5=P(:,5); V6=P(:,6); V7=P(:,7); V8=P(:,8); V9=P(:,9); V10=P(:,10);
h1=P(:,11); h2=P(:,12); h5=P(:,13); h6=P(:,14); h9=P(:,15); h10=P(:,16); 
q=P(:,17); v=P(:,18);

%% Outputs
Vhalf=-30; Vth=-50;
k1=3; k2=8;   
Lth_Ia=60.007; Lth_II=58.457;  
rhov=0.6;
kvIa=6.2; kdIa=2;
kEMGIa=0.06; kEMGII=0.06;
kIb=1; kdII=1.5;
CIa=0.026; CII=0;
a1=60; a2=7;
Femax=-37.7; Fth=3.393;
Lopt=68;

f = zeros(length(tspan),10);  %neuron outputs
for i=1:8
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k1)).^(-1);
end

Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); 
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); le = Le/Lopt;
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
vmf = v.*hf; vme = -v.*he;

Fl = @(l) exp(-(abs((l^2.3-1)/1.26))^1.62); %force dependence on muscle length
Fv = @(vm,l) (vm<0)*(-0.69-0.17*vm)/(vm-0.69)+(vm>=0)*(0.18-(-5.34*l^2+8.41*l-4.7)*vm)/(vm+0.18); %force dependence on velocity
Fp = @(l) 3.5*log(exp((l-1.4)/0.05)+1)-0.02*(exp(-18.7*(l-0.79))-1);  %passive force
Fe = zeros(length(tspan),1);
for i = 1:length(tspan)
    Fe(i) = Femax*(f(i,10).*Fl(le(i)).*Fv(vme(i),le(i))+Fp(le(i))); %extensor muscle force
end

%feedback current 
Iaf = s_Iaf*max(sign(vmf).*kvIa.*abs(vmf/Lth_Ia).^rhov+kdIa*max((Lf-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,9)+CIa,0); 
Iae = s_Iae*max(sign(vme).*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,10)+CIa,0);
Ibe = -s_Ibe*kIb*max(-Fe-Fth,0)/Femax;
IIf = s_IIf*(kdII*max((Lf-Lth_II)/Lth_II,0)+kEMGII*f(:,9)+CII);

%% Plot
figure(1)
xlim=[0 tF];

subplot(5,2,1)
plot(tspan,f(:,1),'-b','LineWidth',1); hold on
plot(tspan,f(:,2),'-r','LineWidth',1); 
ylim_f12=[0,0.573];
shade_st(T0,T0_st,ylim_f12)
axis([xlim ylim_f12]);
ylabel('RG'); set(gca,'FontSize',12)

subplot(5,2,2)
plot(tspan,f(:,5),'-b','LineWidth',1); hold on
plot(tspan,f(:,6),'-r','LineWidth',1); 
ylim_f56=[0,0.5];
shade_st(T0,T0_st,ylim_f56)
axis([xlim ylim_f56]); 
ylabel('PF'); set(gca,'FontSize',12)

subplot(5,2,3)
plot(tspan,f(:,3),'-b','LineWidth',1); hold on
plot(tspan,f(:,4),'-r','LineWidth',1); 
ylim_f34=[0,0.66];
shade_st(T0,T0_st,ylim_f34)
axis([xlim ylim_f34]);
ylabel('In'); set(gca,'FontSize',12)

subplot(5,2,4)
plot(tspan,f(:,9),'-b','LineWidth',1); hold on
plot(tspan,f(:,10),'-r','LineWidth',1); 
ylim_f910=[0,0.42];
shade_st(T0,T0_st,ylim_f910)
axis([xlim ylim_f910]); 
ylabel('Mn'); set(gca,'FontSize',12)

subplot(5,2,5)
plot(tspan,f(:,7),'-b','LineWidth',1); hold on
plot(tspan,f(:,8),'-r','LineWidth',1); 
ylim_f78=[0,0.4852];
shade_st(T0,T0_st,ylim_f78)
axis([xlim ylim_f78]); 
ylabel('Int (blue), Inab-E (red)'); set(gca,'FontSize',12)

subplot(5,2,6)
plot(tspan,q,'-k','LineWidth',1); hold on 
ylim_q=[1.29,1.84];
shade_st(T0,T0_st,ylim_q)
axis([xlim,ylim_q]); 
ylabel('limb angle'); set(gca,'FontSize',12)

subplot(5,2,7)
plot(tspan,Iaf,'-b','LineWidth',1); hold on
ylim_Ia=[0,0.142];
shade_st(T0,T0_st,ylim_Ia)
axis([xlim,ylim_Ia]); 
ylabel('Ia-F'); set(gca,'FontSize',12)

subplot(5,2,8)
plot(tspan,Iae,'-r','LineWidth',1); hold on
ylim_Ia=[0,0.147];
shade_st(T0,T0_st,ylim_Ia)
axis([xlim,ylim_Ia]); 
ylabel('Ia-E'); set(gca,'FontSize',12)

subplot(5,2,9)
plot(tspan,IIf,'-b','LineWidth',1); hold on
ylim_II=[0,0.114];
shade_st(T0,T0_st,ylim_II)
axis([xlim,ylim_II]); 
xlabel('time'); ylabel('II-F'); set(gca,'FontSize',12)

subplot(5,2,10)
plot(tspan,Ibe,'-r','LineWidth',1); hold on
ylim_Ib=[0,0.32];
shade_st(T0,T0_st,ylim_Ib)
axis([xlim,ylim_Ib]); 
xlabel('time'); ylabel('Ib-E'); set(gca,'FontSize',12)