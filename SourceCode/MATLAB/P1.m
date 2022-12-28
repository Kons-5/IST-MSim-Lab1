%% P1
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% Variable definition
K12 = 0.3*3600; K21 = 0.2455*3600; K10 = 0.0643*3600;
V = 3100;                       %V1 = V2
delta = 1000;
N = 3;                          %where N is the number of dosages tested
h = 1;                          %step size
c50 = 7.1903;                   %half maximal effect
T = 21;                         %administration period

%% Init
% Dosage vector  
d = zeros(1,N) + 3;
d = upsample(d,T);

%time vector
size = length(d); t = 1:h:size;

%% define ODE function

% c(t)' juhu= f(c,d)
f = @(c,d)[(1/V*(-K12-K10)*c(1,:) + 1/V*K21*c(2,:) + delta/V*d);
             (1/V*K12*c(1,:) - 1/V*K21*c(2,:))];

%% Function Call

%Euler's integration and plot
Euler(f,t,h,d,size);

%% Function Def
function Euler(f,t,h,d,size)
    
    % Preallocation
    c = zeros(2,size); 

    % Simple Euler's integration
    for i = 1:(size-1)
        c(:,i+1) = c(:,i) + h*f(c(:,i),d(i));
    end
 
    % Plot pharmacocinetics
    figure();

    % Concentration in compartment 1 and 2
    yyaxis left
    plot(t,c(1,:), LineWidth=1.5);
    hold on; plot(t,c(2,:), LineStyle='-', LineWidth=1.5, Color='#D95319'); 
    ylabel('\textbf{Concentra\c{c}\~ao [mg/Kg]}', 'Interpreter', 'latex','FontSize',12);
    
    % Dose administration
    yyaxis right
    stem(t,d,'.', LineWidth=0.25, MarkerSize=14,Color="k"); ylim([-0.6, 3.6]); 
    ylabel('\textbf{Dose Administrada [mg/Kg]}', 'Interpreter', 'latex','FontSize',12, Color='k');
    
    grid, grid minor; xlim([1, 63]);
    xlabel('\textbf{Tempo [dia]}', 'Interpreter','latex','FontSize',12); 
    legend({'\textbf{Concentra\c{c}\~ao no comp. 1}','\textbf{Concentra\c{c}\~ao no comp. 2}'}, 'Interpreter', 'latex', 'Location','south');
    title('\textbf{Resposta do Modelo Farmacocin\''etico (PK)}', 'Interpreter','latex','FontSize',16);

    % change left and right axis to black
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    hold off

end