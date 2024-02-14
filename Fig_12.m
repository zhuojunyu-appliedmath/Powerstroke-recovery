%Fig. 12

clear; clc; close all;

figure

load('FB_strength_corrected.mat');
h1=plot(PS_Iaf(:,2),PS_Iaf(:,1),'.','MarkerSize',14); hold on
h2=plot(PS_IIf(:,2),PS_IIf(:,1),'.','MarkerSize',14);
h3=plot(PS_Iae(:,2),PS_Iae(:,1),'.','MarkerSize',14);
h4=plot(PS_Ibe(:,2),PS_Ibe(:,1),'.','MarkerSize',14);
plot(PS_Iaf(19,2),PS_Iaf(19,1),'.k','MarkerSize',32)
legd=legend([h1 h2 h3 h4],'$\uppercase\expandafter{\romannumeral1}a-F$','$\uppercase\expandafter{\romannumeral2}-F$','$\uppercase\expandafter{\romannumeral1}a-E$','$\uppercase\expandafter{\romannumeral1}b-E$');
set(legd,'Interpreter','latex');
xlabel('Sensitivity'); ylabel('Performance');
axis([0 2.81 .073 .24]);
set(gca,'FontSize',13);
