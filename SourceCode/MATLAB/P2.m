%% P2
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% Function calls
Hill();

%% Function Def
function Hill()
    C50 = 7.1903;
    %create plethera of concentration values
    c2 = 0:0.01:1000;

    %create Hill's equation: Pharmacodynamic model
    u = c2./(C50 + c2);

    %plot pharmacodynamic's function
    semilogx(c2,u), grid; ylim([-0.2, 1.2]);
    
    %assymptotes
    yline(1, '--', '$u_{max}$', 'interpreter','latex', 'LineWidth',1,'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','top',FontSize=12);

    yline(0, '--', 'LineWidth',1,'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','top',FontSize=12);
    
    %labels
    ylabel('Efeito, u','interpreter','latex','FontSize',12); 
    xlabel('Concentra\c{c}\~ao no comp. 2, [mg/kg]','interpreter','latex','FontSize',12);

    title('\textbf{Resposta Modelo F\''armacodin\^amico, [mg/kg]}','interpreter','latex','FontSize',12);
    subtitle('Efeito do F\''armaco em fun\c{c}\~ao em escala semi-logar\''itmica','interpreter','latex','FontSize',12);
   
    
end