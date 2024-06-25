% Fig. 3

clear; clc; close all;

%% Parameters
a0=0.165;               
g=2;                    
F0=10;              
b=4000;  
Ethresh=15;     %CPG synaptic threshold    
gsyn=0.005;     %CPG synaptic conductance   
Efb=-80;        %feedback synaptic threshold
gfb=0.001;      %feedback synaptic conductance
F_ell=2;        %applied load
kappa=1;        %applied load strength
pert=0.1;       %perturbation to the strength
kappa_pert=kappa+pert;  %perturbed strengthdt=0.01; 
dt=0.01;

L0=10; Lslope=1; 

init = [15.0000   19.8248    0.3010    0.7832    0.0000    0.5349    2.6749    0];

%% Find unperturbed solution with kappa and perturbed solution with kappa+pert
[T0,T0_ps,init_ps_u,end_ps_u,init_re_u,end_re_u] = phases(gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,init);  %unperturbed
T0_re=T0-T0_ps;
[Tp,Tp_ps,init_ps_p,end_ps_p,init_re_p,end_re_p] = phases(gsyn,Ethresh,gfb,Efb,kappa+pert,F_ell,L0,Lslope,init);  %perturbed
Tp_re=Tp-Tp_ps;

%% Estimate linear shift in time T1 by calculating lTRC
% for power stroke
dx_in_ps=(init_ps_p'-init_ps_u')/pert;
dx_out_ps=(end_ps_p'-end_ps_u')/pert;
[T1_ps,n_in_ps] = lTRC(gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,dt,T0_ps,init_ps_u,dx_in_ps,dx_out_ps,1,[]);
% for recovery
dx_in_re=(init_re_p'-init_re_u')/pert;
dx_out_re=(end_re_p'-end_re_u')/pert;
[T1_re,~] = lTRC(gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,dt,T0-T0_ps,init_re_u,dx_in_re,dx_out_re,0,n_in_ps);

nu1_ps=T1_ps/T0_ps;   %nu1 in power stroke
nu1_re=T1_re/T0_re;   %nu1 in recovery

T1=T1_ps+T1_re;    %linear shift in total period
beta0 = T0_ps/T0; 
beta1=(T1_ps*T0-T0_ps*T1)/T0^2;  %linear shift in the proportion of power-stroke duration

%% Find the iSRC with piecewise uniform rescaling; start from power stroke
tF=4*T0; tspan=0:dt:tF;
init_src = [init_ps_u 0 (init_ps_p-init_ps_u)/pert];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[Tsrc,G] = ode15s(@iSRC,tspan,init_src,options,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,nu1_ps,nu1_re);
V1=G(:,1); V2=G(:,2); N1=G(:,3); N2=G(:,4); A1=G(:,5); A2=G(:,6); x=G(:,7); y=G(:,8);  %trajectory
src=G(:,9:15);  %isrc

%% Compare unperturbed and perturbed
kappa1=1;
init_ps1=[15.0000   19.8248    0.3010    0.7832    0.0000    0.5349    2.6749   0];
T01=3054.62; T0_ps1=1544.16;
tspan1=0:.01:4*T01;

kappa2=2;
init_ps2=[15.0000   19.7314    0.3006    0.7845   -0.0000    0.5295    3.4569   0];
T02=2830.76; T0_ps2=1366.33;
tspan2=0:.01:T02;

[~,P1] = ode15s(@model,tspan1,init_ps1,[],gsyn,Ethresh,gfb,Efb,kappa1,F_ell,L0,Lslope);
[~,P2] = ode15s(@model,tspan2,init_ps2,[],gsyn,Ethresh,gfb,Efb,kappa2,F_ell,L0,Lslope);


%% Plot trajectory comparison and iSRC

figure
xlim=[0 T01];
xlim_src=[2*T01 3*T01+8];
ind=T01/T02;

subplot(4,2,1)
plot(tspan1,P1(:,1),'-b','LineWidth',1); hold on
plot(tspan1,P1(:,2),'-r','LineWidth',1);
plot(tspan2*ind,P2(:,1),'--b','LineWidth',1); 
plot(tspan2*ind,P2(:,2),'--r','LineWidth',1);
plot(tspan1,Ethresh*ones(1,length(tspan1)),'--m','LineWidth',1.5);
ylim_V=[-48,58];
shade_ps(T01,T0_ps1,ylim_V)
plot(T0_ps2*ind*ones(1,length(ylim_V(1):ylim_V(2))),ylim_V(1):ylim_V(2),':g','LineWidth',1.5);
axis([xlim ylim_V]); 
ylabel('V'); set(gca,'FontSize',12);

subplot(4,2,3)
plot(Tsrc,src(:,1),'-b','LineWidth',1.5); hold on
plot(Tsrc,src(:,2),'-r','LineWidth',1.5);
ylim_srcV=[-11 16];
shade_ps(T0,T0_ps,ylim_srcV)
axis([xlim_src ylim_srcV]); hold off
ylabel('iSRC - V'); set(gca,'FontSize',12);

subplot(4,2,2)
plot(tspan1,P1(:,3),'-b','LineWidth',1); hold on
plot(tspan1,P1(:,4),'-r','LineWidth',1); 
plot(tspan2*ind,P2(:,3),'--b','LineWidth',1); 
plot(tspan2*ind,P2(:,4),'--r','LineWidth',1);
ylim_N=[0.25,0.85];
shade_ps(T01,T0_ps1,ylim_N);
plot(T0_ps2*ind*ones(1,length(ylim_N(1):.01:ylim_N(2))),ylim_N(1):.01:ylim_N(2),':g','LineWidth',1.5);
axis([xlim ylim_N]); 
ylabel('N'); set(gca,'FontSize',12);

subplot(4,2,4)
plot(Tsrc,src(:,3),'-b','LineWidth',1.5); hold on
plot(Tsrc,src(:,4),'-r','LineWidth',1.5); 
ylim_srcN=[-0.05 0.07];
shade_ps(T0,T0_ps,ylim_srcN)
axis([xlim_src ylim_srcN]); hold off
ylabel('iSRC - N'); set(gca,'FontSize',12);

subplot(4,2,5)
plot(tspan1,P1(:,5),'-b','LineWidth',1); hold on
plot(tspan1,P1(:,6),'-r','LineWidth',1);
plot(tspan2*ind,P2(:,5),'--b','LineWidth',1); 
plot(tspan2*ind,P2(:,6),'--r','LineWidth',1);
ylim_A=[-0.09,1.08];
shade_ps(T01,T0_ps1,ylim_A);
plot(T0_ps2*ind*ones(1,length(ylim_A(1):.01:ylim_A(2))),ylim_A(1):.01:ylim_A(2),':g','LineWidth',1.5);
axis([xlim ylim_A]); 
ylabel('A'); set(gca,'FontSize',12);

subplot(4,2,7)
plot(Tsrc,src(:,5),'-b','LineWidth',1.5); hold on
plot(Tsrc,src(:,6),'-r','LineWidth',1.5); 
ylim_srcA=[-.8 .8];
shade_ps(T0,T0_ps,ylim_srcA)
axis([xlim_src ylim_srcA]); hold off
xlabel('time'); ylabel('iSRC - A'); set(gca,'FontSize',12);

subplot(4,2,6)
plot(tspan1,P1(:,7),'-k','LineWidth',1); hold on
plot(tspan2*ind,P2(:,7),'--K','LineWidth',1);
ylim_x=[-1.4,4];
shade_ps(T01,T0_ps1,ylim_x);
plot(T0_ps2*ind*ones(1,length(ylim_x(1):.1:ylim_x(2))),ylim_x(1):.1:ylim_x(2),':g','LineWidth',1.5);
axis([xlim,ylim_x]); 
ylabel('x'); set(gca,'FontSize',12);

subplot(4,2,8)
plot(Tsrc,src(:,7),'-k','LineWidth',1.5); hold on
ylim_srcx=[0.6 1.52];
shade_ps(T0,T0_ps,ylim_srcx)
axis([xlim_src ylim_srcx]); hold off
xlabel('time'); ylabel('iSRC - x'); set(gca,'FontSize',12);