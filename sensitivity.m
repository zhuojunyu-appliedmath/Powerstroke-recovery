% Compute average and sensitivity of average using iSRC

close all;

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
pert=0.05;       %perturbation to the strength
kappa_pert=kappa+pert;  %perturbed strength
dt=0.01; 

Lc=10; Sc=0.5;
Lslope_range = 1; 

init_ps_u = [15.0000000000134,19.8250085179671,0.300957072526219,0.783148580449431,-3.20213847588544e-42,0.534953604170953,2.67423237468234];

Q=zeros(length(Lslope_range),1);
D=zeros(length(Lslope_range),1);
data=zeros(length(Lslope_range),9);

for m = 1:length(Lslope_range)

Lslope=Lslope_range(m);
L0=Lc-Lslope*log(Sc/(1-Sc))/2;

init = [init_ps_u 0];

%% Find unperturbed solution with kappa and perturbed solution with kappa+pert
[T0,T0_ps,init_ps_u,end_ps_u,init_re_u,end_re_u] = phases(gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,init);  %unperturbed
T0_re=T0-T0_ps;
[Tp,Tp_ps,init_ps_p,end_ps_p,init_re_p,end_re_p] = phases(gsyn,Ethresh,gfb,Efb,kappa+pert,F_ell,L0,Lslope,init);  %perturbed
Tp_re=Tp-Tp_ps;

beta0 = T0_ps/T0; %beta_p=Tp_ps/Tp;
%beta1 = (beta_p-beta0)/pert;  %numerical estimate

%% Estimate linear shift in time T1 by calculating lTRC
%T1_ps = (Tp_ps-T0_ps)/pert; % numerical estimate
%T1_re = (Tp_re-T0_re)/pert;  % numerical estimate
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

T1=T1_ps+T1_re;  %linear shift in total period
beta1=(T1_ps*T0-T0_ps*T1)/T0^2;  %linear shift in the proportion of power-stroke duration

%% Find the iSRC with piecewise uniform rescaling; start from power stroke
tF=T0; tspan=0:dt:tF;
init_src = [init_ps_u 0 (init_ps_p-init_ps_u)/pert];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,G] = ode15s(@iSRC,tspan,init_src,options,gsyn,Ethresh,gfb,Efb,kappa,F_ell,L0,Lslope,nu1_ps,nu1_re);
V1=G(:,1); V2=G(:,2); N1=G(:,3); N2=G(:,4); A1=G(:,5); A2=G(:,6); x=G(:,7); y=G(:,8);  %trajectory
src=G(:,9:15);  %isrc

%% Calculate average and sensitivity of the average
tspan_ps=0:dt:T0_ps;
% int_dQdkappa = zeros(length(tspan_ps),1);
% 
% for i = 1:length(tspan_ps)
%     
%     L1 = 10+x(i); L2 = 10-x(i);
%     u1 = (1/2)*V1(i); u2 = (1/2)*V2(i);
%     U1 = (1.03-4.31*exp(-0.198*u1))*(u1>=8);
%     U2 = (1.03-4.31*exp(-0.198*u2))*(u2>=8);
%     LT1 = -3*sqrt(3)*(L1-1)*(L1-5)*(L1-15)/2/625;
%     LT2 = -3*sqrt(3)*(L2-1)*(L2-5)*(L2-15)/2/625;
%     dLT1 = -3*sqrt(3)*(3*L1^2-42*L1+95)/2/625;
%     dLT2 = -3*sqrt(3)*(3*L2^2-42*L2+95)/2/625;
%     a1 = g*(A1(i)-a0); a2 = g*(A2(i)-a0);
%     F1 = (F0*a1*LT1)*(u1>=8)*(a1>=0);
%     F2 = (F0*a2*LT2)*(u2>=8)*(a2>=0);
%          
%     if u1>=8 && a1>=0
%         dF1dA1 = F0*g*LT1;
%         dF1dx = F0*a1*dLT1*0.8;
%     else 
%         dF1dA1 = 0;
%         dF1dx = 0;
%     end
% 
%     if u2>=8 && a2>=0
%         dF2dA2 = F0*g*LT2;
%         dF2dx = F0*a2*dLT2*(-0.8);
%     else 
%         dF2dA2 = 0;
%         dF2dx = 0;
%     end
%     
%     int_Q=(F2-F1+kappa*F_ell)/b;
%     Jq = [0  0  0  0  -dF1dA1/b  dF2dA2/b  (dF2dx-dF1dx)/b];
%     int_dQdkappa(i)=-beta0*Jq*src(i,:)'-beta0*F_ell/b-beta1*int_Q;
%     
% end

delta_x =  y(floor(T0_ps/dt));
delta_x1 = src(floor(T0_ps/dt),7)-src(1,7);
D(m) = abs(delta_x1/delta_x-T1/T0);

Q(m) = delta_x/T0;
%dQdkappa_1 = abs(dt*trapz(int_dQdkappa)/T0_ps);
dQdkappa_2 = Q(m)*D(m);

data(m,:)=[Q(m) T0 T1 T1/T0 delta_x delta_x1 delta_x1/delta_x D(m) dQdkappa_2];

end

init_ps_u
init_re_u

%% Plot solution
figure(1) %trajectory
xlim=[0 tF];

subplot(3,2,1)
plot(tspan,V1,'-b','LineWidth',1); hold on
plot(tspan,V2,'-r','LineWidth',1);
plot(tspan,Ethresh*ones(1,length(tspan)),'--m','LineWidth',1);
ylim_V=[-48,55];
shade_ps(T0,T0_ps,ylim_V)
axis([xlim ylim_V]); %hold off
xlabel('time'); ylabel('V'); set(gca,'FontSize',11);

subplot(3,2,2)
plot(tspan,N1,'-b','LineWidth',1); hold on
plot(tspan,N2,'-r','LineWidth',1); 
ylim_N=[0.28,0.8];
shade_ps(T0,T0_ps,ylim_N)
axis([xlim ylim_N]); %hold off
xlabel('time'); ylabel('N'); set(gca,'FontSize',11);

subplot(3,2,3)
plot(tspan,A1,'-b','LineWidth',1); hold on
plot(tspan,A2,'-r','LineWidth',1)
ylim_A=[-0.09,1.08];
shade_ps(T0,T0_ps,ylim_A)
axis([xlim ylim_A]); %hold off
xlabel('time'); ylabel('A'); set(gca,'FontSize',11);

subplot(3,2,4)
plot(tspan,x,'--k','LineWidth',1); hold on
ylim_x=[-3,3];
shade_ps(T0,T0_ps,ylim_x)
axis([xlim,ylim_x]); %hold off
xlabel('time'); ylabel('x'); set(gca,'FontSize',11);

subplot(3,2,5)
plot(tspan,y,'-k','LineWidth',1); hold on
ylim_y=[-1 8];
shade_ps(T0,T0_ps,ylim_y)
axis([xlim ylim_y]); %hold off
xlabel('time'); ylabel('y'); set(gca,'FontSize',11);

subplot(3,2,6)
plot(V1,N1,'-b','LineWidth',1); hold on
plot(Ethresh*ones(1,length([0.18:0.01:0.77])),[0.18:0.01:0.77],'--m','LineWidth',1); %hold off
axis([-48,55,0.28,0.8])
xlabel('V_1'); ylabel('N_1'); set(gca,'FontSize',11);

figure(2)  %iSRC
subplot(2,2,1)
plot(tspan,src(:,1),'-b','LineWidth',1.5); hold on
plot(tspan,src(:,2),'-r','LineWidth',1.5);
ylim_srcV=[-7 10];
shade_ps(T0,T0_ps,ylim_srcV)
axis([xlim ylim_srcV]); hold off
xlabel('time'); ylabel('iSRC - V'); set(gca,'FontSize',12);

subplot(2,2,2)
plot(tspan,src(:,3),'-b','LineWidth',1.5); hold on
plot(tspan,src(:,4),'-r','LineWidth',1.5); 
ylim_srcN=[-0.04 0.05];
shade_ps(T0,T0_ps,ylim_srcN)
axis([xlim ylim_srcN]); hold off
xlabel('time'); ylabel('iSRC - N'); set(gca,'FontSize',12);

subplot(2,2,3)
plot(tspan,src(:,5),'-b','LineWidth',1.5); hold on
plot(tspan,src(:,6),'-r','LineWidth',1.5); 
ylim_srcA=[-.5 .5];
shade_ps(T0,T0_ps,ylim_srcA)
axis([xlim ylim_srcA]); hold off
xlabel('time'); ylabel('iSRC - A'); set(gca,'FontSize',12);

subplot(2,2,4)
plot(tspan,src(:,7),'-k','LineWidth',1.5); hold on
ylim_srcx=[0.45 1.52];
shade_ps(T0,T0_ps,ylim_srcx)
axis([xlim ylim_srcx]); hold off
xlabel('time'); ylabel('iSRC - x'); set(gca,'FontSize',12);
