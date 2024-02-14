% Fig. 8

clear; clc; close all;

figure

load('excitation-contralateral.mat');
plot(dec_9(41:end,9),dec_9(41:end,1),'.b','MarkerSize',15); hold on
plot(dec_10(41:end,9),dec_10(41:end,1),'.k','MarkerSize',15);
plot(dec_11(41:end,9),dec_11(41:end,1),'.r','MarkerSize',15);
plot(inc_9(41:end,9),inc_9(41:end,1),'+b','MarkerSize',7.5);
plot(inc_10(41:end,9),inc_10(41:end,1),'+k','MarkerSize',7.5);
plot(inc_11(41:end,9),inc_11(41:end,1),'+r','MarkerSize',7.5);
plot(inc_11(end,9),inc_11(end,1),'.g','MarkerSize',28);

load('inhibition-contralateral.mat');
plot(dec_9(40:end,9),dec_9(40:end,1),'.b','MarkerSize',15); hold on
plot(dec_10(40:end,9),dec_10(40:end,1),'.k','MarkerSize',15);
plot(dec_11(40:end,9),dec_11(40:end,1),'.r','MarkerSize',15);
plot(inc_9(40:end,9),inc_9(40:end,1),'+b','MarkerSize',7.5);
plot(inc_10(40:end,9),inc_10(40:end,1),'+k','MarkerSize',7.5);
plot(inc_11(40:end,9),inc_11(40:end,1),'+r','MarkerSize',7.5);
plot(inc_11(end,9),inc_11(end,1),'.g','MarkerSize',28);

xlabel('Sensitivity'); ylabel('Performance');
title('Contralateral Feedback: top - excitation, bottom - inhibition');
axis([.00015 .00027 .001196 .001375]);
set(gca,'FontSize',13);