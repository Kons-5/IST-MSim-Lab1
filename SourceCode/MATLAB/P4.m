%% P3
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% Function calls

% init dosage vector 
T = 21; N = 24;
d = zeros(1,N) + 3; %<- bolus therapy dosage
d = upsample(d,T);   

% init concentration
c = concentration(d);

% init effect
u = Hill(c);

% init volume growth
v = Tumor(u);

% plot
PlotRelevantVariables(c, d, v, u);

% find optimal administration period
Teste = Optimal();

%% Function defs
function c = concentration(d)

    K12 = 0.3*3600; K21 = 0.2455*3600; K10 = 0.0643*3600;
    V = 3100;                       %V1 = V2
    delta = 1000;
    h = 1;                          %step size           
    
    %create f function -> c(t)' = f(c,d)
    f = @(c,d)[(1/V*(-K12-K10)*c(1,:) + 1/V*K21*c(2,:) + delta/V*d);
                 (1/V*K12*c(1,:) - 1/V*K21*c(2,:))];
    
    %preallocation
    size = length(d);
    c = zeros(2,size); 

    %simple Euler's integration
    for i = 1:(size-1)
        c(:,i+1) = c(:,i) + h*f(c(:,i),d(i));
    end
end
%% 
function v = Tumor(u)

    %tumor growth variables
    a = 0.09; b = 1; Kt=10; h = 1; p0 = 1; %1mm^3 
    size = length(u); 

    %create logistic equation
    l = @(v,u)a*v*(1-v./Kt) - b*u.*v;

    %preallocating for speed
    v = zeros(1,size); 
    v(1) = p0;

    %simple Euler's integration
    for i = 1:(size - 1)
        v(i+1) = v(i) + h*l(v(i),u(i));
    end
end
%%
function u = Hill(c)

    %half maximal effect
    c50 = 7.1903;
    %create Hill's equation: Pharmacodynamic model
    u = c(2,:)./(c50 + c(2,:));
end
%%
function PlotRelevantVariables(c, d, v, u)

    %time vector
    size = length(d); t = 1:size;

    figure()
    subplot(2,1,1);
    title('\textbf{Simula\c{c}\~ao do Sistema completo}','FontSize',16,'Interpreter', 'latex');

    % Concentration in compartment 1 and 2
    yyaxis left
    plot(t,c(1,:), LineWidth=1.5);
    hold on; plot(t,c(2,:), LineStyle='-', LineWidth=1.5, Color='#D95319'); 
    ylabel('\textbf{Concentra\c{c}\~ao, [mg/kg]}', 'Interpreter', 'latex','FontSize',12);

    % Dose administration
    yyaxis right
    stem(t,d,'.',Color='k', LineWidth=0.5); ylim([-0.6, 3.6]); 
    ylabel('\textbf{Dose Administrada, [mg/Kg]}', 'Interpreter', 'latex','FontSize',12, Color='k');

    %change left and right axis to black
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend({'\textbf{c1 - Concentra\c{c}\~ao no comp. 1}','\textbf{c2 - Concentra\c{c}\~ao no comp. 2}', '\textbf{Dose administrada}'}, ...
        'Interpreter', 'latex', 'Location','south');
    xlim([1,300]); grid, grid minor;
    
    subplot(2,1,2)

    %tumor volume progression plot
    yyaxis left
    plot(t,v, LineWidth=1.5); ylabel('\textbf{Volume, [mm$^3$]}','interpreter','latex','FontSize',12); ylim([0, 2.5])
    yyaxis right; plot(t,u, LineWidth=1.5); ylabel('\textbf{Efeito, u}','interpreter','latex','FontSize',12); ylim([0, 0.25])

    xlabel('\textbf{Tempo [dias]}','interpreter','latex','FontSize',12)
    xlim([1,300]); grid, grid minor;

    legend({'\textbf{Volume}','\textbf{Efeito}'}, ...
        'Interpreter', 'latex', 'Location','south');
end
%%
function T = Optimal()

    opVolume = 0.10; N = 25; d = zeros(1,N) + 3;

    for T = 1:21

    clear dosage; dosage = upsample(d,T); 

    % Calculate volume evolution
    c = concentration(dosage);
    u = Hill(c);
    v = Tumor(u);

    % evaluate volume at 25 day mark
    if v(25) >= opVolume
        break
    end
    end 
    fprintf("Optimal periodicity: %i days\n", T);
    Plot();
end
%%
function Plot()
    N = 25; d = zeros(1,N) + 3;

    for T = [2, 3, 4]

    clear dosage; dosage = upsample(d,T); 
    size = length(dosage); t = 1:size;

    % Calculate volume evolution
    c = concentration(dosage);
    u = Hill(c);
    v = Tumor(u);
    if T == 3
        s=plot(t,v,'LineStyle','-',LineWidth=1.5, Color='b');
        plot(25,v(25),'.',MarkerSize=22, Color='b');
        datatip(s,25,0.12)
        continue
    end
    figure(2);
    plot(t,v,'LineStyle','--',LineWidth=1.5);
    hold on
    end
    figure(2);
    xlim([1,25]);
    xlabel('\textbf{Tempo [dias]}','interpreter','latex','FontSize',12)
    ylabel('\textbf{Volume, [mm$^3$]}','interpreter','latex','FontSize',12);

    legend({'\textbf{$T = 2$}','\textbf{$T = 3$}','','\textbf{$T = 4$}'}, ...
        'Interpreter', 'latex', 'Location','best');

    title('\textbf{Volume do tumor}','\textbf{consoante o espa\c{c}amento entre administra\c{c}\~oes}','FontSize',16,'Interpreter', 'latex');
end