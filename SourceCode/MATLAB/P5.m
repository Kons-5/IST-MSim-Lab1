%% P5
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;

%% Function calls
N = 12; T = 21;

% init dosage vector -> up front dosing
d = dVariablePeriod(N);
c = concentration(d);
u = Hill(c);
v = Tumor(u);
PlotRelevantVariables(d, v);

% init dosage vector -> metronomic pattern
d = dMetronomicPattern();
c = concentration(d);
u = Hill(c);
v = Tumor(u);
PlotRelevantVariables(d, v);

% init dosage vector -> metronomic amplitude comparison
comparisonMetronome(T)


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
function u = Hill(c)

    %half maximal effect
    c50 = 7.1903;
    %create Hill's equation: Pharmacodynamic model
    u = c(2,:)./(c50 + c(2,:));
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
function d = dVariablePeriod(N)
    d = [];    
    pos = 1;
    d(1) = 3;
    for ii = 1:N
        for jj = 1:fibonacci(ii)
            pos = pos + 1;
            d(pos) = 0;
        end
    
        pos = pos + 1; %apply dose after spacing
        d(pos) = 3;
    end
end
%%
function comparisonMetronome(T)

    d = zeros(1,40) + 3; % -> make amplitude constant
    d = upsample(d, T+1);

    c = concentration(d);
    u = Hill(c);
    v = Tumor(u);

    %time vector
    size = length(d); t = 1:size;
    
    %tumor volume progression plot
    figure();
    plot(t,v, LineWidth=1.5);
    hold on

    % Dose administration
    d = d./2.15;
    stem(t,d,'.', LineWidth=0.25, MarkerSize=14,Color="#0072BD");  
    hold on

    d = repmat([6 0 3],1,13);
    d = upsample(d, T+1);

    %time vector
    size = length(d); t = 1:size;

    c = concentration(d);
    u = Hill(c);
    v = Tumor(u);

    %tumor volume progression plot
    plot(t,v, LineWidth=1.5); 
    hold on

    % Dose administration
    d = d./2.15; grid; grid minor;
    stem(t,d,'.', LineWidth=0.25, MarkerSize=14,Color="#EDB120");  
    ylabel('\textbf{Volume, [mm$^3$]}', 'Interpreter', 'latex','FontSize',12);
    xlabel('\textbf{Tempo, [dias]}', 'Interpreter', 'latex','FontSize',12)
    title('\textbf{Amplitude Vari\''avel vs Amplitude Constante}','Interpreter', 'latex','FontSize',16)
    hold off

    xlim([0, 340]);ylim([0,3]);

end
%%
function d = dMetronomicPattern()
    d = repmat([3 0 0 0 3 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0],1,31);
end
%%
function PlotRelevantVariables(d, v)

    %time vector
    size = length(d); t = 1:size;

    figure()
    
    %tumor volume progression plot
    plot(t,v, LineWidth=1.5); ylabel('\textbf{Volume, [mm$^3$]}','interpreter','latex','FontSize',12); ylim([0, 1.5])
    hold on

    % Dose administration
    d = d./2.15;
    stem(t,d,'.',Color='k', LineWidth=0.5, MarkerSize=12);  
    ylabel('\textbf{Volume, [mm$^3$]}', 'Interpreter', 'latex','FontSize',12);
    hold off 

    xlabel('\textbf{Tempo [dias]}','interpreter','latex','FontSize',12)
    xlim([1,244]); grid, grid minor;

    legend({'\textbf{Volume}','\textbf{Dose administrada (3 mg)}'}, ...
        'Interpreter', 'latex', 'Location','south');

    title('\textbf{Metronomic Dosing: Padr\~ao Fibonacci}','FontSize',16,'Interpreter', 'latex');
end
