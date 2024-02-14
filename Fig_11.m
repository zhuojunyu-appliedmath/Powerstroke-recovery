% Fig. 11

% Neuron numbering: 1, RG-F; 2, RG-E; 3, In-F; 4, In-E; 5, PF-F; 6, PF-E; 
% 7, Int; 8, Inab-E; 9, Mn-F; 10, Mn-E; 

clear; clc; close all;

%feedback strength
s_Iaf=1;
s_Iae=1;
s_Ibe=1;
s_IIf=1;

%% Find unperturbed solution with kappa and perturbed solution with kappa+pert
kappa1=0;
init_st1=[-64.8809361683031,-36.9594183760335,-58.6472831890731,-27.7473050183001,-62.2355295469095,-33.0647778390756,-63.6072508638440,-31.9808035411915,-63.9624843671704,-32.3185157661491,0.533521095965132,0.301667837348812,0.392933826237958,0.167281516515830,0.403832356834355,0.232028801749969,1.29920016313509,4.57591682652905e-17];
T01=1035.2913; T0_st1=719.0293; tspan1=0:.1:T01;

kappa2=0.01;
init_st2=[-64.8789926315448,-36.9606806071470,-58.6412380001848,-27.7559112872561,-62.2323595287412,-33.0676796829582,-63.6048161752362,-31.9872839402687,-63.9619273124349,-32.3222652330766,0.531587947976039,0.301620956446874,0.393517298296791,0.166100251202825,0.404832514459309,0.231023701752848,1.29980596200574,4.34506726117775e-18];
T02=1104.3230; T0_st2=791.2377; tspan2=0:.1:T02; 

[~,P1] = ode15s(@model,tspan1,init_st1,[],kappa1,s_Iaf,s_Iae,s_Ibe,s_IIf);
[~,P2] = ode15s(@model,tspan2,init_st2,[],kappa2,s_Iaf,s_Iae,s_Ibe,s_IIf);


%% Outputs
Vhalf=-30; Vth=-50;
k1=3;       %for Mn
k2=8;       %otherwise
ls=300;

f1 = zeros(length(tspan1),10);  %neuron outputs
for i=1:2
    f1(:,i) = (P1(:,i)>=Vth).*(1+exp(-(P1(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f1(:,i) = (P1(:,i)>=Vth).*(1+exp(-(P1(:,i)-Vhalf)/k1)).^(-1);
end

f2 = zeros(length(tspan2),10);  %neuron outputs
for i=1:2
    f2(:,i) = (P2(:,i)>=Vth).*(1+exp(-(P2(:,i)-Vhalf)/k2)).^(-1);
end
for i=9:10
    f2(:,i) = (P2(:,i)>=Vth).*(1+exp(-(P2(:,i)-Vhalf)/k1)).^(-1);
end

y1 = zeros(length(tspan1),1);
for i = 1:floor(T0_st1/0.1)
    y1(i) = -ls*(cos(P1(i,17)-kappa1)-cos(P1(1,17)-kappa1));
end
y1(floor(T0_st1/0.1)+1:end)=y1(floor(T0_st1/0.1));
Per1=y1(end)/T01;

y2 = zeros(length(tspan2),1);
for i = 1:floor(T0_st2/0.1)
    y2(i) = -ls*(cos(P2(i,17)-kappa1)-cos(P2(1,17)-kappa1));
end
y2(floor(T0_st2/0.1)+1:end)=y2(floor(T0_st2/0.1));
Per2=y2(end)/T02;

%% Plot
figure(1)
xlim=[0 T01];
ind=T01/T02;

subplot(2,2,1)
plot(tspan1,f1(:,1),'-b','LineWidth',1); hold on
plot(tspan1,f1(:,2),'-r','LineWidth',1); 
plot(tspan2*ind,f2(:,1),'--b','LineWidth',1); 
plot(tspan2*ind,f2(:,2),'--r','LineWidth',1); 
ylim_f12=[0,0.573];
plot(T0_st2*ind*ones(1,length(ylim_f12(1):.01:ylim_f12(2))),ylim_f12(1):.01:ylim_f12(2),':g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_f12)
axis([xlim ylim_f12]);
ylabel('RG output'); set(gca,'FontSize',12)

subplot(2,2,3)
plot(tspan1,f1(:,9),'-b','LineWidth',1); hold on
plot(tspan1,f1(:,10),'-r','LineWidth',1); 
plot(tspan2*ind,f2(:,9),'--b','LineWidth',1);
plot(tspan2*ind,f2(:,10),'--r','LineWidth',1); 
ylim_f910=[0,0.42];
plot(T0_st2*ind*ones(1,length(ylim_f910(1):.01:ylim_f910(2))),ylim_f12(1):.01:ylim_f910(2),':g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_f910)
axis([xlim ylim_f910]); 
xlabel('time'); ylabel('Mn output'); set(gca,'FontSize',12)

subplot(2,2,2)
plot(tspan1,P1(:,17),'-k','LineWidth',1); hold on 
plot(tspan2*ind,P2(:,17),'--k','LineWidth',1);
ylim_q=[1.29,1.84];
plot(T0_st2*ind*ones(1,length(ylim_q(1):.01:ylim_q(2))),ylim_q(1):.01:ylim_q(2),':g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_q)
axis([xlim,ylim_q]); 
ylabel('limb angle'); set(gca,'FontSize',12)

subplot(2,2,4)
plot(tspan1,y1,'-k','LineWidth',1); hold on 
plot(tspan2*ind,y2,'--k','LineWidth',1);
ylim_y=[-10,165];
plot(T0_st2*ind*ones(1,length(ylim_y(1):ylim_y(2))),ylim_y(1):ylim_y(2),':g','LineWidth',1.5);
shade_st(T01,T0_st1,ylim_y)
axis([xlim,ylim_y]); 
xlabel('time'); ylabel('progress'); set(gca,'FontSize',12)
