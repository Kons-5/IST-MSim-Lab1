%% P4
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% Function calls
% exemplify area of mutation
AreaOfMutation();

%concentration for therapy in use
cp = BolusTherapy();

%calculate resistance vector
resistance = SimulateResistance(cp);

%exemplify competitive resistance aquisition through double bolus therapy
[PharmaDyR,PharmadyNR, RC50] = CompetitiveHill(resistance, cp);
volumeNResistance = Tumor(PharmadyNR);
volumeComp = Tumor(PharmaDyR);

%plot competitive
PlotGraphsC(cp, PharmaDyR, PharmadyNR, RC50);

%exemplify non competitive resistance aquisition """"
[PharmaDyR, RC50] = NCompetitiveHill(resistance, cp);
volumeNComp = Tumor(PharmaDyR);

%plot noncompetitive
PlotGraphsNC(cp, PharmaDyR, PharmadyNR);

%plot tumor Growth
TumorGrowthPlot(volumeNResistance,volumeComp,volumeNComp);



%% Function def
function AreaOfMutation()


    %polynomial representing drug concentration
    d = @(x)2.64*x.^4-23.92*x.^3+77.28*x.^2-104*x+45.4;
    e = -2:0.005:5;
   

    %plot area where mutation occurs
    below = d(e); 
    below = below.*(below<=0);
    figure();
    area(e,below, 'LineStyle', 'none','FaceAlpha',.3);
    hold on;

    plot(e,d(e), 'LineWidth',2,'Color','b');
    ylim([-5,5]); xlim([0,4.5]);

    %treshold line
    yline(0, '-', '$L_r$', 'interpreter','latex', 'LineWidth',2,'FontSize',15);
    
    set(gca, 'YTick', []);
    xticks(0:4);
    xlabel('\textbf{Tempo, [dias]}', 'interpreter','latex','FontSize',12)
    ylabel('\textbf{Concentra\c{c}\~ao de Atezolizumab, [mg/Kg]}', 'interpreter','latex', 'FontSize',12)
    title('\textbf{Exemplo do Periodo de Aquisi\c{c}\~ao de Resist\^encia}', 'interpreter','latex', 'FontSize',16)
end
%%
function cp = BolusTherapy()

    K12 = 0.3*3600; K21 = 0.2455*3600; K10 = 0.0643*3600;
    V = 3100;                       %V1 = V2
    delta = 1000;
    N = 24;                         %where N is the number of dosages tested
    h = 1;                          %step size
    T = 21;                         %administration period
    
    %init dosage vector  
    d = zeros(1,N) + 20; %<- bolus therapy dosage
    d = upsample(d,T); size = length(d);                    
    
    %create f function -> c(t)' = f(c,d)
    f = @(c,d)[(1/V*(-K12-K10)*c(1,:) + 1/V*K21*c(2,:) + delta/V*d);
                 (1/V*K12*c(1,:) - 1/V*K21*c(2,:))];

    %preallocation
    c = zeros(2,size); 

    %simple Euler's integration
    for i = 1:size
        c(:,i+1) = c(:,i) + h*f(c(:,i),d(i));
    end

    cp = c(2,:);
end
%%
function resistance = SimulateResistance(cp)
    
    Kr = 0.1; Lr = 4.5; h = 1;
    f = @(cp)(Kr*max(0,Lr - cp));
    
    %preallocation
    resistance = zeros(1,length(cp));

    %simple Euler's integration
    for i = 1:(length(cp) - 1)
        resistance(i+1) = resistance(i) + h*f(cp(i));
    end
end
%% 
function PlotGraphsC(cp,PharmaDyR,PharmadyNR, RC50)

    Lr = 4.5; 

    %init time vector
    size = length(cp); t = 1:size;

    figure(); 
    subplot(2,1,1); 
    title('\textbf{Modelo de Resist\^encia (Competitivo)}', 'Interpreter','latex','FontSize',16);

    %concentration in compartment of effect
    yyaxis left
    plot(t,cp, 'LineWidth',1.5); ylim([0, 6.5]); ylabel('\textbf{Concentra\c{c}\~ao no comp. 2, [mg/kg]}', ...
        'FontSize',12, 'LineWidth',1.5,'Interpreter','latex','FontSize',12);

    %treshold line
    yline(Lr, '--', '$L_r$', 'interpreter','latex', 'LineWidth',1,'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','top',FontSize=15);

    %C50 shift right
    yyaxis right
    plot(t,RC50, 'LineWidth',1.5); ylabel('\textbf{C50, [mg/kg]}','FontSize',12,'Interpreter','latex');

    legend({'\textbf{Concentra\c{c}\~ao no comp. 2, [mg/kg]}','\textbf{\textit{Treshold}}', '\textbf{C50, [mg/kg]}'},'Interpreter','latex')

    xlim([0,128]); grid; grid minor;
  
    %change in effect
    subplot(2,1,2)
    yyaxis right; plot(t,PharmaDyR, 'LineWidth',1.5); ylim([0, 0.5]); ylabel('\textbf{Efeito com Resist\^encia}','FontSize',12,'Interpreter','latex');
    yyaxis left; plot(t,PharmadyNR, 'LineWidth',1.5); ylim([0, 0.5]); ylabel('\textbf{Efeito sem Resist\^encia}','FontSize',12,'Interpreter','latex');

    legend({'\textbf{Efeito com Resist\^encia}', '\textbf{Efeito sem Resist\^encia}'},'Interpreter','latex')

    xlim([0,128]); xlabel('\textbf{Tempo, [dias]}','Interpreter','latex','FontSize',12); grid; grid minor;
end
%%
function [PharmadyR,PharmadyNR, RC50] = CompetitiveHill(resistance, cp)

    %shift right with effect
    C50 = 7.1903;
    RC50 = (1 + resistance).*C50;

    %%create Hill's equation
    PharmadyR = cp./(RC50 + cp);
    PharmadyNR = cp./(C50 + cp);    
end
%%
function [PharmaDyR, RC50] = NCompetitiveHill(resistance, cp)    

    %shift right with effect
    RC50 = 7.1903;
    
    %Non competitive resistance model
    f = 1./(1 + resistance);

    %%create Hill's equation
    PharmaDyR = f.*((cp.*f)./(RC50 + cp));
end
%%
function PlotGraphsNC(cp, PharmaDyR, PharmadyNR)
    Lr = 4.5; C50 = 7.1903;

    %init time vector
    size = length(cp); t = 1:size;

    figure();
    subplot(2,1,1);
    title('\textbf{Modelo de Resist\^encia (N\~ao Competitivo)}', 'Interpreter','latex','FontSize',16);

    %concentration in compartment of effect
    yyaxis left
    plot(t,cp, 'LineWidth',1.5); ylim([0, 8]); ylabel('\textbf{Concentra\c{c}\~ao no comp. 2, [mg/kg]}', ...
        'FontSize',12, 'LineWidth',1.5,'Interpreter','latex','FontSize',12);

    %treshold line
    yline(Lr, '--', '$L_r$', 'interpreter','latex', 'LineWidth',1,'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','top',FontSize=15);

    %C50 line
    yyaxis right
    yline(C50, '-', 'LineWidth',1,'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','top',FontSize = 15, LineWidth = 1.5, Color='#D95319'); 
    ylim([0, 8]); ylabel('\textbf{C50, [mg/kg]}','FontSize',12,'Interpreter','latex');

    legend({'\textbf{Concentra\c{c}\~ao no comp. 2, [mg/kg]}','\textbf{\textit{Treshold}}', '\textbf{C50, [mg/kg]}'},'Interpreter','latex')

   xlim([0,128]); grid; grid minor;

    %change in effect
    subplot(2,1,2)
    yyaxis right; plot(t,PharmaDyR, 'LineWidth',1.5); ylim([0, 0.14]); ylabel('\textbf{Efeito com Resist\^encia}','FontSize',12,'Interpreter','latex');
    yyaxis left; plot(t,PharmadyNR, 'LineWidth',1.5); ylim([0, 0.5]); ylabel('\textbf{Efeito sem Resist\^encia}','FontSize',12,'Interpreter','latex'); 

    legend({'\textbf{Efeito com Resist\^encia}', '\textbf{Efeito sem Resist\^encia}'},'Interpreter','latex')

    xlim([0,128]); xlabel('Tempo, [dias]','FontSize',12); grid; grid minor;
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
function TumorGrowthPlot(volumeNResistance, volumeComp, volumeNComp)

    size = length(volumeNResistance); t = 1:size;
    
    figure();
    plot(t,volumeNResistance, 'LineWidth',1.5);
    hold on; plot(t,volumeComp, 'LineWidth',1.5);
    hold on; plot(t,volumeNComp, 'LineWidth',1.5);

    hold off;

    legend({'\textbf{Evolu\c{c}\~ao do Volume sem Resist\^encia}', '\textbf{Evolu\c{c}\~ao do volume com Antagonista Competitivo}', ...
        '\textbf{Evolu\c{c}\~ao do Volume com Antagonista N\~ao Competitivo}'}, 'Interpreter', 'latex','Location','best');

    xlabel('\textbf{Tempo, [dias]}','Interpreter', 'latex','FontSize',12);
    ylabel('\textbf{Volume, [mm$^{3}$]}','Interpreter', 'latex','FontSize', 12);

    title('\textbf{Evolu\c{c}\~ao do crescimento do tumor}','Interpreter', 'latex','FontSize', 16);

    xlim([0, 150]); grid, grid minor;
end