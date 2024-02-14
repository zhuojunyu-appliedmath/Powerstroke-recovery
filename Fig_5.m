% Fig. 5

clear; clc; close all;

%% Parameters
a0=0.165;               
g=2;                    
F0=10;              
b=4000;  
Ethresh=15;     %CPG synaptic threshold    
gsyn=0.005;     %CPG synaptic conductance   
gfb=0.001;      %feedback synaptic conductance
F_ell=2;        %applied load
kappa=1;        %applied load strength
pert=0.05;       %perturbation to the strength
kappa_pert=kappa+pert;  %perturbed strength
dt=10;

L0=10; Lslope=50000; %almost constant feeedback activation

%% Solutions
% excitatory feedback
init_ps_ex=[15.0000   22.0509    0.3787    0.7681   -0.0000    0.6468    2.3759    0];
T0_ex=2488.75; T0_ps_ex=1264.31;
tspan_ex=0:dt:T0_ex;
[~,P_ex] = ode15s(@model,tspan_ex,init_ps_ex,[],gsyn,Ethresh,gfb,80,kappa,F_ell,L0,Lslope);
V1_ex=P_ex(:,1); V2_ex=P_ex(:,2); N1_ex=P_ex(:,3); N2_ex=P_ex(:,4); A1_ex=P_ex(:,5); A2_ex=P_ex(:,6); x_ex=P_ex(:,7); y_ex=P_ex(:,8);

% inhibitory feedback
init_ps_in=[15.0000   19.3518    0.3217    0.7650   -0.0000    0.5065    2.4907    0];
T0_in=2798.38; T0_ps_in=1410.89;
tspan_in=0:dt:T0_in;
[~,P_in] = ode15s(@model,tspan_in,init_ps_in,[],gsyn,Ethresh,gfb,-80,kappa,F_ell,L0,Lslope);
V1_in=P_in(:,1); V2_in=P_in(:,2); N1_in=P_in(:,3); N2_in=P_in(:,4); A1_in=P_in(:,5); A2_in=P_in(:,6); x_in=P_in(:,7); y_in=P_in(:,8);

%% Plot
figure
xlim=[0 T0_in];

subplot(3,2,1)
plot(tspan_in,V1_in,'-b','LineWidth',1); hold on
plot(tspan_in,V2_in,'-r','LineWidth',1);
plot(tspan_ex,V1_ex,'--b','LineWidth',1);
plot(tspan_ex,V2_ex,'--r','LineWidth',1);
plot(tspan_in,Ethresh*ones(1,length(tspan_in)),'--m','LineWidth',1.5);
ylim_V=[-48,58];
shade_ps(T0_in,T0_ps_in,ylim_V);
plot(T0_ps_ex*ones(1,length(ylim_V(1):.1:ylim_V(2))),ylim_V(1):.1:ylim_V(2),':g','LineWidth',1.5);
plot(T0_ex*ones(1,length(ylim_V(1):.1:ylim_V(2))),ylim_V(1):.1:ylim_V(2),'--g','LineWidth',1.5);
axis([xlim ylim_V]); 
xlabel('time'); ylabel('V'); set(gca,'FontSize',12);

subplot(3,2,2)
plot(tspan_in,N1_in,'-b','LineWidth',1); hold on
plot(tspan_in,N2_in,'-r','LineWidth',1); 
plot(tspan_ex,N1_ex,'--b','LineWidth',1); 
plot(tspan_ex,N2_ex,'--r','LineWidth',1); 
ylim_N=[0.25,0.83];
shade_ps(T0_in,T0_ps_in,ylim_N);
plot(T0_ps_ex*ones(1,length(ylim_N(1):.01:ylim_N(2))),ylim_N(1):.01:ylim_N(2),':g','LineWidth',1.5);
plot(T0_ex*ones(1,length(ylim_N(1):.01:ylim_N(2))),ylim_N(1):.01:ylim_N(2),'--g','LineWidth',1.5);
axis([xlim ylim_N]); 
xlabel('time'); ylabel('N'); set(gca,'FontSize',12);

subplot(3,2,3)
plot(tspan_in,A1_in,'-b','LineWidth',1); hold on
plot(tspan_in,A2_in,'-r','LineWidth',1)
plot(tspan_ex,A1_ex,'--b','LineWidth',1); 
plot(tspan_ex,A2_ex,'--r','LineWidth',1)
ylim_A=[-0.09,1.08];
shade_ps(T0_in,T0_ps_in,ylim_A);
plot(T0_ps_ex*ones(1,length(ylim_A(1):.01:ylim_A(2))),ylim_A(1):.01:ylim_A(2),':g','LineWidth',1.5);
plot(T0_ex*ones(1,length(ylim_A(1):.01:ylim_A(2))),ylim_A(1):.01:ylim_A(2),'--g','LineWidth',1.5);
axis([xlim ylim_A]); 
xlabel('time'); ylabel('A'); set(gca,'FontSize',12);

subplot(3,2,4)
plot(tspan_in,x_in,'-k','LineWidth',1); hold on
plot(tspan_ex,x_ex,'--k','LineWidth',1); 
ylim_x=[-1.4,3];
shade_ps(T0_in,T0_ps_in,ylim_x);
plot(T0_ps_ex*ones(1,length(ylim_x(1):.1:ylim_x(2))),ylim_x(1):.1:ylim_x(2),':g','LineWidth',1.5);
plot(T0_ex*ones(1,length(ylim_x(1):.1:ylim_x(2))),ylim_x(1):.1:ylim_x(2),'--g','LineWidth',1.5);
axis([xlim,ylim_x]); 
xlabel('time'); ylabel('x'); set(gca,'FontSize',12);

subplot(3,2,5)
plot(tspan_in,y_in,'-k','LineWidth',1); hold on
plot(tspan_ex,y_ex,'--k','LineWidth',1); 
ylim_y=[-.5,4.1];
shade_ps(T0_in,T0_ps_in,ylim_y);
plot(T0_ps_ex*ones(1,length(ylim_y(1):.1:ylim_y(2))),ylim_y(1):.1:ylim_y(2),':g','LineWidth',1.5);
plot(T0_ex*ones(1,length(ylim_y(1):.1:ylim_y(2))),ylim_y(1):.1:ylim_y(2),'--g','LineWidth',1.5);
axis([xlim ylim_y]);
xlabel('time'); ylabel('Progress'); set(gca,'FontSize',12);

subplot(3,2,6)
plot(V1_in,N1_in,'.b','MarkerSize',5); hold on
plot(V1_ex,N1_ex,'.r','MarkerSize',5); 
plot(Ethresh*ones(1,length([0.25:0.01:0.83])),[0.25:0.01:0.83],'--m','LineWidth',1.5); 
axis([-45,55,0.25,0.83])
xlabel('V_1'); ylabel('N_1'); set(gca,'FontSize',12);
