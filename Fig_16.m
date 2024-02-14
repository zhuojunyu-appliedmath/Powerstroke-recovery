% Fig. 16

clear; clc; close all;

load s_Ibe.mat
figure
h1=plot(s_Ibe(:,1),s_Ibe(:,6),'.b','MarkerSize',12); hold on
h2=plot(s_Ibe(:,1),s_Ibe(:,9),'.r','MarkerSize',12); 
plot(s_Ibe(46,1),s_Ibe(46,6),'.k','MarkerSize',35); 
axis([0.8 7 -20 6]);
legend([h1 h2],'y_1/y_0','T_1/T_0');
xlabel('s_{Ib-E}'); ylabel('Timing or shape change');
set(gca,'FontSize',13);