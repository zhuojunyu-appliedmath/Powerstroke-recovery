% Fig. 15

% Neuron numbering: 1, RG-F; 2, RG-E; 3, In-F; 4, In-E; 5, PF-F; 6, PF-E; 
% 7, Int; 8, Inab-E; 9, Mn-F; 10, Mn-E; 

clear; clc; close all;

%% Parameters
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
w34=0.44; w44=0.484;
ls=300;

Fl = @(l) exp(-(abs((l^2.3-1)/1.26))^1.62); %force dependence on muscle length
Fv = @(vm,l) (vm<0)*(-0.69-0.17*vm)/(vm-0.69)+(vm>=0)*(0.18-(-5.34*l^2+8.41*l-4.7)*vm)/(vm+0.18); %force dependence on velocity
Fp = @(l) 3.5*log(exp((l-1.4)/0.05)+1)-0.02*(exp(-18.7*(l-0.79))-1);  %passive force

s_Iae=1;
s_Ibe=1;
s_IIf=1;

kappa=0;

%% s_Iaf=1
s_Iaf1=1;
init_st1=[-64.8809378117064,-36.9594157966064,-58.6472908751815,-27.7472995924213,-62.2355322784133,-33.0647767211889,-63.6072521709523,-31.9808006359821,-63.9624843683281,-32.3185145652701,0.533520879043655,0.301667961622656,0.392933844976166,0.167281647228994,0.403832124325481,0.232029016617407,1.29920012370126,8.29921889103028e-19];
T01=1035.22; T0_st1=719.03; tspan1=0:.1:T01;
[~,P1] = ode15s(@model,tspan1,init_st1,[],kappa,s_Iaf1,s_Iae,s_Ibe,s_IIf);
q=P1(:,17); v=P1(:,18);

f1 = zeros(length(tspan1),10);  %neuron outputs
for i=3:4
    f1(:,i) = (P1(:,i)>=Vth).*(1+exp(-(P1(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f1(:,i) = (P1(:,i)>=Vth).*(1+exp(-(P1(:,i)-Vhalf)/k1)).^(-1);
end

Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); 
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); le = Le/Lopt;
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
vmf = v.*hf; vme = -v.*he;
Fe = zeros(length(tspan1),1);
for i = 1:length(tspan1)
    Fe(i) = Femax*(f1(i,10).*Fl(le(i)).*Fv(vme(i),le(i))+Fp(le(i))); %extensor muscle force
end

Iae = s_Iae*max(sign(vme).*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f1(:,10)+CIa,0);
Ibe = -s_Ibe*kIb*max(-Fe-Fth,0)/Femax;

FB1 = w34*Iae+w44*Ibe;

%% s_Iaf=1.1
s_Iaf2=1.1;
init_st2=[-64.8777844733178,-36.8471367545377,-58.3972436036739,-27.7301317447472,-62.2330736136349,-33.0374318949528,-63.6122338608154,-31.9800212319549,-63.9623783672757,-32.2990959429003,0.521764490336965,0.308060261652990,0.382620397609172,0.171954323686490,0.393029428140901,0.240124262560802,1.31300541372862,1.32215437379733e-16];
T02=964.88; T0_st2=657.92; tspan2=0:.1:T02;
[~,P2] = ode15s(@model,tspan2,init_st2,[],kappa,s_Iaf2,s_Iae,s_Ibe,s_IIf);
q=P2(:,17); v=P2(:,18);

f2 = zeros(length(tspan2),10);  %neuron outputs
for i=3:4
    f2(:,i) = (P2(:,i)>=Vth).*(1+exp(-(P2(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f2(:,i) = (P2(:,i)>=Vth).*(1+exp(-(P2(:,i)-Vhalf)/k1)).^(-1);
end

Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); 
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); le = Le/Lopt;
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
vmf = v.*hf; vme = -v.*he;
Fe = zeros(length(tspan2),1);
for i = 1:length(tspan2)
    Fe(i) = Femax*(f2(i,10).*Fl(le(i)).*Fv(vme(i),le(i))+Fp(le(i))); %extensor muscle force
end

Iae = s_Iae*max(sign(vme).*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f2(:,10)+CIa,0);
Ibe = -s_Ibe*kIb*max(-Fe-Fth,0)/Femax;

FB2 = w34*Iae+w44*Ibe;

%% Performance
y1 = zeros(length(tspan1),1);
for i = 1:floor(T0_st1/0.1)
    y1(i) = -ls*(cos(P1(i,17)-kappa)-cos(P1(1,17)-kappa));
end
y1(floor(T0_st1/0.1)+1:end)=y1(floor(T0_st1/0.1));
Per1=y1(end)/T01;

y2 = zeros(length(tspan2),1);
for i = 1:floor(T0_st2/0.1)
    y2(i) = -ls*(cos(P2(i,17)-kappa)-cos(P2(1,17)-kappa));
end
y2(floor(T0_st2/0.1)+1:end)=y2(floor(T0_st2/0.1));
Per2=y2(end)/T02;

%% Plot
figure(1)
xlim=[0 T01];

subplot(2,2,1)
plot(tspan1,f1(:,3),'-b','LineWidth',1); hold on
plot(tspan1,f1(:,4),'-r','LineWidth',1); 
plot(tspan2,f2(:,3),'--b','LineWidth',1); 
plot(tspan2,f2(:,4),'--r','LineWidth',1); 
ylim_f34=[0,0.66];
plot(T0_st2*ones(1,length(ylim_f34(1):.01:ylim_f34(2))),ylim_f34(1):.01:ylim_f34(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_f34(1):.01:ylim_f34(2))),ylim_f34(1):.01:ylim_f34(2),'--g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_f34)
axis([xlim ylim_f34]);
ylabel('In output'); set(gca,'FontSize',12)

subplot(2,2,2)
plot(tspan1,P1(:,17),'-k','LineWidth',1); hold on 
plot(tspan2,P2(:,17),'--k','LineWidth',1);
ylim_q=[1.29,1.84];
plot(T0_st2*ones(1,length(ylim_q(1):.01:ylim_q(2))),ylim_q(1):.01:ylim_q(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_q(1):.01:ylim_q(2))),ylim_q(1):.01:ylim_q(2),'--g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_q)
axis([xlim,ylim_q]); 
ylabel('limb angle'); set(gca,'FontSize',12)

subplot(2,2,3)
plot(tspan1,FB1,'-r','LineWidth',1); hold on
plot(tspan2,FB2,'--r','LineWidth',1); 
ylim_ef=[0,0.23];
plot(T0_st2*ones(1,length(ylim_ef(1):.01:ylim_ef(2))),ylim_ef(1):.01:ylim_ef(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_ef(1):.01:ylim_ef(2))),ylim_ef(1):.01:ylim_ef(2),'--g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_ef)
axis([xlim,ylim_ef]); 
xlabel('time'); ylabel('FB to In-E'); set(gca,'FontSize',12)

subplot(2,2,4)
plot(tspan1,y1,'-k','LineWidth',1); hold on 
plot(tspan2,y2,'--k','LineWidth',1);
ylim_y=[-10,165];
plot(T0_st2*ones(1,length(ylim_y(1):ylim_y(2))),ylim_y(1):ylim_y(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_y(1):ylim_y(2))),ylim_y(1):ylim_y(2),'--g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_y)
axis([xlim,ylim_y]); 
xlabel('time'); ylabel('progress'); set(gca,'FontSize',12)
