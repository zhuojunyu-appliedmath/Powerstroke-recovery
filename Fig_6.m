% Fig. 6

clear; clc; close all;

%% Parameters
Ethresh=15;     %CPG synaptic threshold    
gsyn=0.005;     %CPG synaptic conductance   
Efb=-80;        %feedback synaptic threshold (inhibition -80; excitation +80)
gfb=0.001;      %feedback synaptic conductance
F_ell=2;        %applied load
kappa=1;        %applied load strength
dt=0.01;
L0=11;        

%% Solve ODE system
Lslope1=0.5; 
init_ps1=[15.0000   20.3119    0.3005    0.7773   -0.0000    0.5623    2.8590   0];
T01=3027.50; T0_ps1=1547.17; tspan1=0:dt:T01;
[~,P1] = ode15s(@model,tspan1,init_ps1,[],gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope1);
V11=P1(:,1); V21=P1(:,2); N11=P1(:,3); N21=P1(:,4); A11=P1(:,5); A21=P1(:,6); x1=P1(:,7); y1=P1(:,8);

Lslope2=2; 
init_ps2=[15.0000   19.7953    0.3020    0.7760   -0.0000    0.5332    2.7025   0];
T02=2996.98; T0_ps2=1510.63; tspan2=0:dt:T02;
[~,P2] = ode15s(@model,tspan2,init_ps2,[],gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope2);
V12=P2(:,1); V22=P2(:,2); N12=P2(:,3); N22=P2(:,4); A12=P2(:,5); A22=P2(:,6); x2=P2(:,7); y2=P2(:,8);

%% Muscle forces
a0=0.165;               
g=2;                    
F0=10;              
b=4000;   

LT11=zeros(length(tspan1),1); LT21=zeros(length(tspan1),1);
F11=zeros(length(tspan1),1); F21=zeros(length(tspan1),1);

for i=1:length(tspan1)
    
V1=V11(i); V2=V21(i); A1=A11(i); A2=A21(i); x=x1(i);

L1 = 10+x; L2 = 10-x;
LT11(i) = -3*sqrt(3)*(L1-1)*(L1-5)*(L1-15)/2/625;
LT21(i) = -3*sqrt(3)*(L2-1)*(L2-5)*(L2-15)/2/625;
a1 = g*(A1-a0); a2 = g*(A2-a0);
u1 = (1/2)*V1; u2 = (1/2)*V2;
F11(i) = (F0*a1*LT11(i))*(u1>=8)*(a1>=0);
F21(i) = (F0*a2*LT21(i))*(u2>=8)*(a2>=0);

end

LT12=zeros(length(tspan2),1); LT22=zeros(length(tspan2),1);
F12=zeros(length(tspan2),1); F22=zeros(length(tspan2),1);

for i=1:length(tspan2)
    
V1=V12(i); V2=V22(i); A1=A12(i); A2=A22(i); x=x2(i);

L1 = 10+x; L2 = 10-x;
LT12(i) = -3*sqrt(3)*(L1-1)*(L1-5)*(L1-15)/2/625;
LT22(i) = -3*sqrt(3)*(L2-1)*(L2-5)*(L2-15)/2/625;
a1 = g*(A1-a0); a2 = g*(A2-a0);
u1 = (1/2)*V1; u2 = (1/2)*V2;
F12(i) = (F0*a1*LT12(i))*(u1>=8)*(a1>=0);
F22(i) = (F0*a2*LT22(i))*(u2>=8)*(a2>=0);

end

%% Plot trajectory
figure(1)
xlim=[0 T01];

subplot(4,2,1)
plot(tspan1,V11,'-b','LineWidth',1); hold on
plot(tspan1,V21,'-r','LineWidth',1);
plot(tspan2,V12,'--b','LineWidth',1); 
plot(tspan2,V22,'--r','LineWidth',1);
plot(tspan1,Ethresh*ones(1,length(tspan1)),'--m','LineWidth',1.5);
ylim_V=[-45,53];
shade_ps(T01,T0_ps1,ylim_V);
plot(T0_ps2*ones(1,length(ylim_V(1):ylim_V(2))),ylim_V(1):ylim_V(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_V(1):ylim_V(2))),ylim_V(1):ylim_V(2),'--g','LineWidth',1.5);
axis([xlim ylim_V]); 
xlabel('time'); ylabel('V'); set(gca,'FontSize',12);

subplot(4,2,2)
plot(tspan1,N11,'-b','LineWidth',1); hold on
plot(tspan1,N21,'-r','LineWidth',1); 
plot(tspan2,N12,'--b','LineWidth',1); 
plot(tspan2,N22,'--r','LineWidth',1); 
ylim_N=[0.25,0.82];
shade_ps(T01,T0_ps1,ylim_N);
plot(T0_ps2*ones(1,length(ylim_N(1):.01:ylim_N(2))),ylim_N(1):.01:ylim_N(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_N(1):.01:ylim_N(2))),ylim_N(1):.01:ylim_N(2),'--g','LineWidth',1.5);
axis([xlim ylim_N]); 
xlabel('time'); ylabel('N'); set(gca,'FontSize',12);

subplot(4,2,3)
plot(tspan1,A11,'-b','LineWidth',1); hold on
plot(tspan1,A21,'-r','LineWidth',1);
plot(tspan2,A12,'--b','LineWidth',1); 
plot(tspan2,A22,'--r','LineWidth',1);
ylim_A=[-0.09,1.08];
shade_ps(T01,T0_ps1,ylim_A);
plot(T0_ps2*ones(1,length(ylim_A(1):.01:ylim_A(2))),ylim_A(1):.01:ylim_A(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_A(1):.01:ylim_A(2))),ylim_A(1):.01:ylim_A(2),'--g','LineWidth',1.5);
axis([xlim ylim_A]); 
xlabel('time'); ylabel('A'); set(gca,'FontSize',12);

subplot(4,2,4)
plot(tspan1,x1,'-k','LineWidth',1); hold on
plot(tspan2,x2,'--k','LineWidth',1); hold on
ylim_x=[-1.1,3.1];
shade_ps(T01,T0_ps1,ylim_x);
plot(T0_ps2*ones(1,length(ylim_x(1):.1:ylim_x(2))),ylim_x(1):.1:ylim_x(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_x(1):.1:ylim_x(2))),ylim_x(1):.1:ylim_x(2),'--g','LineWidth',1.5);
axis([xlim,ylim_x]); 
xlabel('time'); ylabel('x'); set(gca,'FontSize',12);

subplot(4,2,5)
plot(tspan1,LT11,'-b','LineWidth',1); hold on
plot(tspan1,LT21,'-r','LineWidth',1); 
plot(tspan2,LT12,'--b','LineWidth',1); 
plot(tspan2,LT22,'--r','LineWidth',1); 
ylim_LT=[0.38 1.05];
shade_ps(T01,T0_ps1,ylim_LT);
axis([xlim,ylim_LT]); 
xlabel('time'); ylabel('Length-tension'); set(gca,'FontSize',12);

subplot(4,2,6)
plot(tspan1,-F11,'-b','LineWidth',1); hold on
plot(tspan1,F21,'-r','LineWidth',1); 
plot(tspan2,-F12,'--b','LineWidth',1); 
plot(tspan2,F22,'--r','LineWidth',1); 
ylim_F=[-22 22];
shade_ps(T01,T0_ps1,ylim_F);
axis([xlim,ylim_F]); 
xlabel('time'); ylabel('Muscle forces'); set(gca,'FontSize',12);

subplot(4,2,7)
L=7:.1:13;
S1=0.5*(1-tanh((L-L0)/Lslope1)); 
S2=0.5*(1-tanh((L-L0)/Lslope2)); 
plot(L,S1,'-k','LineWidth',1); hold on
plot(L,S2,'--k','LineWidth',1);
plot(L0*ones(length(-.1:.01:.5),1),-.1:.01:.5,':m','LineWidth',1.5);
axis([7 13 -.1 1.1]);
xlabel('muscle length L'); ylabel('S_\infty^{FB}(L)'); set(gca,'FontSize',12);

subplot(4,2,8)
plot(tspan1,y1,'-k','LineWidth',1); hold on
plot(tspan2,y2,'--k','LineWidth',1);
ylim_y=[-.2,3.8];
shade_ps(T01,T0_ps1,ylim_y);
plot(T0_ps2*ones(1,length(ylim_y(1):.1:ylim_y(2))),ylim_y(1):.1:ylim_y(2),':g','LineWidth',1.5);
plot(T02*ones(1,length(ylim_y(1):.1:ylim_y(2))),ylim_y(1):.1:ylim_y(2),'--g','LineWidth',1.5);
axis([xlim ylim_y]);
xlabel('time'); ylabel('Progress'); set(gca,'FontSize',12);