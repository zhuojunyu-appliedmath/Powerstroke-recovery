% Fig. 2

clear; clc; close all;

%% Parameters
Ethresh=15;     %CPG synaptic threshold    
gsyn=0.005;     %CPG synaptic conductance   
Efb=-80;        %feedback synaptic threshold (inhibition -80; excitation +80)
gfb=0.001;      %feedback synaptic conductance
F_ell=2;        %applied load
kappa=1;        %applied load strength

Lslope=1; L0=10;  %sigmoid feedback parameters

init_ps=[15.0000   19.8248    0.3010    0.7832    0.0000    0.5349    2.6749   0];
T0=3054.62; T0_ps=1544.16;

dt=0.01; tF=2*T0; tspan=0:dt:tF;

%% Solve ODE system
[~,P] = ode15s(@model,tspan,init_ps,[],gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope);
V1=P(:,1); V2=P(:,2); N1=P(:,3); N2=P(:,4); A1=P(:,5); A2=P(:,6); x=P(:,7); y=P(:,8);

%% Plot trajectory
figure(1)
xlim=[0 tF];

subplot(3,2,1)
plot(tspan,V1,'-b','LineWidth',1); hold on
plot(tspan,V2,'-r','LineWidth',1);
plot(tspan,Ethresh*ones(1,length(tspan)),'--m','LineWidth',1.5);
ylim_V=[-48,58];
shade_ps(T0,T0_ps,ylim_V)
axis([xlim ylim_V]); 
xlabel('time'); ylabel('V'); set(gca,'FontSize',12);

subplot(3,2,2)
plot(tspan,N1,'-b','LineWidth',1); hold on
plot(tspan,N2,'-r','LineWidth',1); 
ylim_N=[0.25,0.85];
shade_ps(T0,T0_ps,ylim_N)
axis([xlim ylim_N]); 
xlabel('time'); ylabel('N'); set(gca,'FontSize',12);

subplot(3,2,3)
plot(tspan,A1,'-b','LineWidth',1); hold on
plot(tspan,A2,'-r','LineWidth',1)
ylim_A=[-0.09,1.08];
shade_ps(T0,T0_ps,ylim_A)
axis([xlim ylim_A]); 
xlabel('time'); ylabel('A'); set(gca,'FontSize',12);

subplot(3,2,4)
plot(tspan,x,'-k','LineWidth',1); hold on
ylim_x=[-1.4,3.1];
shade_ps(T0,T0_ps,ylim_x)
axis([xlim,ylim_x]); 
xlabel('time'); ylabel('x'); set(gca,'FontSize',12);

subplot(3,2,5)
L=7:.1:13;
S=0.5*(1-tanh((L-L0)/Lslope)); 
plot(L,S,'-k','LineWidth',1); hold on
plot(L0*ones(length(-.1:.01:.5),1),-.1:.01:.5,':m','LineWidth',1.5);
axis([7 13 -.1 1.1]);
xlabel('muscle length L'); ylabel('S_\infty^{FB}(L)'); set(gca,'FontSize',12);

subplot(3,2,6)
plot(tspan,y,'-k','LineWidth',1); hold on
ylim_y=[-.8,8.5];
shade_ps(T0,T0_ps,ylim_y)
axis([xlim ylim_y]);
xlabel('time'); ylabel('Progress'); set(gca,'FontSize',12);