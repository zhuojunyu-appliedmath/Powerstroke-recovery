% Compute average and sensitivity of average using central difference

clear; clc; close all;

%% Parameters
ls=300;
kappa=0;          %unperturbed slope
pert1=0.005;      %perturbation to the slope
pert2=-pert1;     %perturbation to the slope
kappa_pert1=kappa+pert1;  %perturbed slope
kappa_pert2=kappa+pert2;  %perturbed slope
%feedback strength
s_Iaf=1;
s_Iae=1;
%s_Ibe=1;
s_IIf=1;

s_Ibe_range=6.8;

dt=0.1;  %time step
init_st_u=[-66.4736614758742,-34.9464843268594,-58.9164889306804,-13.8775455443294,-64.6784118709875,-29.6038696679355,-65.5885916605917,-18.1576619491843,-63.9653782318916,-28.6551121057005,0.571994275332856,0.242935938989262,0.454100832645719,0.111285978331943,0.386153603606666,0.139525051297081,1.26070113703791,1.33496693826463e-18];

%Per = zeros(length(s_Ibe_range),1); Sen = zeros(length(s_Ibe_range),1);
DATA = zeros(length(s_Ibe_range),9);

for k = 1:length(s_Ibe_range)
    
s_Ibe = s_Ibe_range(k);

%% Find unperturbed solution with kappa and perturbed solution with kappa_pert
[T0,T0_st,init_st_u,~,~,~] = phases(init_st_u,kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);   %unperturbed
[Tp1,Tp_st1,init_st_p1,~,~,~] = phases(init_st_u,kappa_pert1,s_Iaf,s_Iae,s_Ibe,s_IIf); %perturbed
[Tp2,Tp_st2,init_st_p2,~,~,~] = phases(init_st_u,kappa_pert2,s_Iaf,s_Iae,s_Ibe,s_IIf); %perturbed

tspan_u=0:dt:T0_st; 
[~,Pu] = ode15s(@model,tspan_u,init_st_u,[],kappa,s_Iaf,s_Iae,s_Ibe,s_IIf);
q_u=Pu(:,17); 

tspan_p1=0:dt:Tp_st1;
[~,Pp1] = ode15s(@model,tspan_p1,init_st_p1,[],kappa_pert1,s_Iaf,s_Iae,s_Ibe,s_IIf);
q_p1=Pp1(:,17); 

tspan_p2=0:dt:Tp_st2;
[~,Pp2] = ode15s(@model,tspan_p2,init_st_p2,[],kappa_pert2,s_Iaf,s_Iae,s_Ibe,s_IIf);
q_p2=Pp2(:,17); 

%% Performance and sensitivity
y0=ls*(cos(q_u(1))-cos(q_u(end)));
yp1=ls*(cos(q_p1(1)-kappa_pert1)-cos(q_p1(end)-kappa_pert1));
yp2=ls*(cos(q_p2(1)-kappa_pert2)-cos(q_p2(end)-kappa_pert2));
Qu=y0/T0;
Qp1=yp1/Tp1;
Qp2=yp2/Tp2;

y1=(yp1-yp2)/(pert1-pert2);
T1=(Tp1-Tp2)/(pert1-pert2);
dQu=(Qp1-Qp2)/(pert1-pert2);  %central difference

%Per(k) = Qu;   %performance
%Sen(k) = abs(dQu);  %sensitivity

DATA(k,:) = [Qu  abs(dQu)  y0  y1  y1/y0  T0  T1  T1/T0  abs(Qu*(y1/y0-T1/T0))];

end