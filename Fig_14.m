% Fig. 14

clear; clc; close all

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
w13=0.27; w23=0.1566;
w34=0.44; w44=0.484;
ls=300;

Fl = @(l) exp(-(abs((l^2.3-1)/1.26))^1.62); %force dependence on muscle length
Fv = @(vm,l) (vm<0)*(-0.69-0.17*vm)/(vm-0.69)+(vm>=0)*(0.18-(-5.34*l^2+8.41*l-4.7)*vm)/(vm+0.18); %force dependence on velocity
Fp = @(l) 3.5*log(exp((l-1.4)/0.05)+1)-0.02*(exp(-18.7*(l-0.79))-1);  %passive force

s_Iae=1;
s_Ibe=1;
s_IIf=1;

kappa=0;

%% s_Iaf=0.6
s_Iaf=0.6;

% Solution
init_st=[-64.9338394202533,-37.5032917570111,-59.2929904629905,-27.5849548230146,-62.3145236648714,-33.1666516182790,-63.6519965686033,-31.8287053872033,-63.9582134493339,-32.3614717539528,0.595161179952789,0.270346834235027,0.452697517956963,0.152828438006310,0.461616032800943,0.194902527064849,1.18973432447685,7.53680652044016e-19];
T0=1921.77; T0_st=1551.54;
[tspan,P] = ode15s(@model,[0:.1:T0],init_st,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
V1=P(:,1); V2=P(:,2); V3=P(:,3); V4=P(:,4); V9=P(:,9); V10=P(:,10); 
h1=P(:,11); h2=P(:,12); 
q=P(:,17); v=P(:,18);

% Output
f = zeros(length(tspan),10);  %neuron outputs
for i=1:8
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k1)).^(-1);
end

% Feedback

Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); 
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); le = Le/Lopt;
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
vmf = v.*hf; vme = -v.*he;
Fe = zeros(length(tspan),1);
for i = 1:length(tspan)
    Fe(i) = Femax*(f(i,10).*Fl(le(i)).*Fv(vme(i),le(i))+Fp(le(i))); %extensor muscle force
end

Iaf = s_Iaf*max(sign(vmf).*kvIa.*abs(vmf/Lth_Ia).^rhov+kdIa*max((Lf-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,9)+CIa,0); 
Iae = s_Iae*max(sign(vme).*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,10)+CIa,0);
Ibe = -s_Ibe*kIb*max(-Fe-Fth,0)/Femax;
IIf = s_IIf*(kdII*max((Lf-Lth_II)/Lth_II,0)+kEMGII*f(:,9)+CII);

% Plot figure

figure(1)
xlim=[0 T0];

subplot(4,2,1)
plot(tspan,f(:,1),'-b','LineWidth',1); hold on
plot(tspan,f(:,2),'-r','LineWidth',1); 
ylim_f12=[0,0.6];
shade_st(T0,T0_st,ylim_f12)
axis([xlim ylim_f12]);
ylabel('RG'); set(gca,'FontSize',12)

subplot(4,2,2)
plot(tspan,f(:,3),'-b','LineWidth',1); hold on
plot(tspan,f(:,4),'-r','LineWidth',1); 
ylim_f34=[0,0.7];
shade_st(T0,T0_st,ylim_f34)
axis([xlim ylim_f34]);
ylabel('In'); set(gca,'FontSize',12)

subplot(4,2,3)
plot(tspan,Iaf,'-b','LineWidth',1); hold on
ylim_Iaf=[0,0.122];
shade_st(T0,T0_st,ylim_Iaf)
axis([xlim,ylim_Iaf]); 
ylabel('Ia-F'); set(gca,'FontSize',12)

subplot(4,2,4)
plot(tspan,Iae,'-r','LineWidth',1); hold on
ylim_Iae=[0,0.175];
shade_st(T0,T0_st,ylim_Iae)
axis([xlim,ylim_Iae]); 
ylabel('Ia-E'); set(gca,'FontSize',12)

subplot(4,2,5)
plot(tspan,IIf,'-b','LineWidth',1); hold on
ylim_II=[0,0.162];
shade_st(T0,T0_st,ylim_II)
axis([xlim,ylim_II]); 
ylabel('II-F'); set(gca,'FontSize',12)

subplot(4,2,6)
plot(tspan,Ibe,'-r','LineWidth',1); hold on
ylim_Ib=[0,0.38];
shade_st(T0,T0_st,ylim_Ib)
axis([xlim,ylim_Ib]); 
ylabel('Ib-E'); set(gca,'FontSize',12)

subplot(4,2,7)
plot(tspan,w13*Iaf+w23*IIf,'-b','LineWidth',1); hold on
ylim_ff=[0,0.058];
shade_st(T0,T0_st,ylim_ff)
axis([xlim,ylim_ff]); 
xlabel('time'); ylabel('FB to In-F'); set(gca,'FontSize',12)

subplot(4,2,8)
plot(tspan,w34*Iae+w44*Ibe,'-r','LineWidth',1); hold on
ylim_ef=[0,0.26];
shade_st(T0,T0_st,ylim_ef)
axis([xlim,ylim_ef]); 
xlabel('time'); ylabel('FB to In-E'); set(gca,'FontSize',12)


%% s_Iaf=0.63
s_Iaf=0.63;

% Solution
init_st=[-64.9321588967424,-37.4230544541875,-59.2526956990025,-27.5894939780623,-62.3114731426295,-33.1543186885797,-63.6508223541506,-31.8417127452162,-63.9592014424777,-32.3545743489699,0.586286053848836,0.274865674251350,0.442659714187110,0.153003577608644,0.451615195034297,0.199769144123560,1.20161882840168,4.13418385058001e-17];
T0=1716.66; T0_st=1352.91;
[tspan,P] = ode15s(@model,[0:.1:T0],init_st,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
V1=P(:,1); V2=P(:,2); V3=P(:,3); V4=P(:,4); V9=P(:,9); V10=P(:,10); 
h1=P(:,11); h2=P(:,12); 
q=P(:,17); v=P(:,18);

% Output
f = zeros(length(tspan),10);  %neuron outputs
for i=1:8
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f(:,i) = (P(:,i)>=Vth).*(1+exp(-(P(:,i)-Vhalf)/k1)).^(-1);
end

% Feedback

Lf = sqrt(a1^2+a2^2-2*a1*a2*cos(q)); 
Le = sqrt(a1^2+a2^2-2*a1*a2*cos(pi-q)); le = Le/Lopt;
hf = a1*a2*sin(q)./Lf; he = a1*a2*sin(pi-q)./Le;      
vmf = v.*hf; vme = -v.*he;
Fe = zeros(length(tspan),1);
for i = 1:length(tspan)
    Fe(i) = Femax*(f(i,10).*Fl(le(i)).*Fv(vme(i),le(i))+Fp(le(i))); %extensor muscle force
end

Iaf = s_Iaf*max(sign(vmf).*kvIa.*abs(vmf/Lth_Ia).^rhov+kdIa*max((Lf-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,9)+CIa,0); 
Iae = s_Iae*max(sign(vme).*kvIa.*abs(vme/Lth_Ia).^rhov+kdIa*max((Le-Lth_Ia)/Lth_Ia,0)+kEMGIa*f(:,10)+CIa,0);
Ibe = -s_Ibe*kIb*max(-Fe-Fth,0)/Femax;
IIf = s_IIf*(kdII*max((Lf-Lth_II)/Lth_II,0)+kEMGII*f(:,9)+CII);

% Plot figure

figure(1)

subplot(4,2,1)
plot(tspan,f(:,1),'--b','LineWidth',1); 
plot(tspan,f(:,2),'--r','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_f12(1):.05:ylim_f12(2))),ylim_f12(1):.05:ylim_f12(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_f12(1):.05:ylim_f12(2))),ylim_f12(1):.05:ylim_f12(2),'--g','LineWidth',1.5);

subplot(4,2,2)
plot(tspan,f(:,3),'--b','LineWidth',1);
plot(tspan,f(:,4),'--r','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_f34(1):.05:ylim_f34(2))),ylim_f34(1):.05:ylim_f34(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_f34(1):.05:ylim_f34(2))),ylim_f34(1):.05:ylim_f34(2),'--g','LineWidth',1.5);

subplot(4,2,3)
plot(tspan,Iaf,'--b','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_Iaf(1):.01:ylim_Iaf(2))),ylim_Iaf(1):.01:ylim_Iaf(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_Iaf(1):.01:ylim_Iaf(2))),ylim_Iaf(1):.01:ylim_Iaf(2),'--g','LineWidth',1.5);

subplot(4,2,4)
plot(tspan,Iae,'--r','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_Iae(1):.01:ylim_Iae(2))),ylim_Iae(1):.01:ylim_Iae(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_Iae(1):.01:ylim_Iae(2))),ylim_Iae(1):.01:ylim_Iae(2),'--g','LineWidth',1.5);

subplot(4,2,5)
plot(tspan,IIf,'--b','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_II(1):.01:ylim_II(2))),ylim_Iae(1):.01:ylim_II(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_II(1):.01:ylim_II(2))),ylim_Iae(1):.01:ylim_II(2),'--g','LineWidth',1.5);

subplot(4,2,6)
plot(tspan,Ibe,'--r','LineWidth',1);
plot(T0_st*ones(1,length(ylim_Ib(1):.01:ylim_Ib(2))),ylim_Ib(1):.01:ylim_Ib(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_Ib(1):.01:ylim_Ib(2))),ylim_Ib(1):.01:ylim_Ib(2),'--g','LineWidth',1.5);

subplot(4,2,7)
plot(tspan,w13*Iaf+w23*IIf,'--b','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_ff(1):.002:ylim_ff(2))),ylim_ff(1):.002:ylim_ff(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_ff(1):.002:ylim_ff(2))),ylim_ff(1):.002:ylim_ff(2),'--g','LineWidth',1.5);

subplot(4,2,8)
plot(tspan,w34*Iae+w44*Ibe,'--r','LineWidth',1); 
plot(T0_st*ones(1,length(ylim_ef(1):.01:ylim_ef(2))),ylim_ef(1):.01:ylim_ef(2),':g','LineWidth',1.5);
plot(T0*ones(1,length(ylim_ef(1):.01:ylim_ef(2))),ylim_ef(1):.01:ylim_ef(2),'--g','LineWidth',1.5);
