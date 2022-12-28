%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% y = k + ae^(bt)
a = 1;
b = log(5/a + 1)/(23);
k = -a;
t = 0:0.5:35;
y = k + a*exp(b*t);

%% plot shenanigans
p = plot(t,y, 'LineWidth', 2), grid, grid minor; hold on;
xlim([0 35]); ylim([0 15]); yticks([0 1 2 3 4 5 10 15]);

xline(23, '-.'); yline(5, '-.');
plot(23, 5, 'r.', 'MarkerSize', 32);
datatip(p, 23, 5, 'Location', 'southeast');

%labels
xlabel('\textbf{Concentra\c{c}\~ao do f\''armaco, [mg/kg]}','interpreter','latex','FontSize',12);
ylabel('\textbf{N\''ivel de toxicidade}','interpreter','latex','FontSize',12);   
title('\textbf{Curva de toxicidade}', 'Interpreter','latex','FontSize',16);