% Fig. 7
clear; clc; close all;

figure

load('excitation-contralateral.mat');
h1=plot(log10(Lslope(21:end)),dec_9(21:end,4),'-r','LineWidth',2); hold on
h2=plot(log10(Lslope(21:end)),dec_9(21:end,7),'--r','LineWidth',2);

load('inhibition-contralateral.mat');
h3=plot(log10(Lslope_9(21:end)),dec_9(21:end,4),'-b','LineWidth',2); 
h4=plot(log10(Lslope_9(21:end)),dec_9(21:end,7),'--b','LineWidth',2); 

legend([h4 h2 h1 h3],'Inhibitory, y_1/y_0','Excitatory, y_1/y_0','Excitatory, T_1/T_0','Inhibitory, y_1/y_0');
axis([-0.22 4.7 -.1 .33]);
xlabel('log(L_{slope})'); ylabel('Timing or shape change');
set(gca,'FontSize',13);
