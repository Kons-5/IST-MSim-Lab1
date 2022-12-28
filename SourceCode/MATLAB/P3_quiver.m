%% P3: Quiver
%% clear->erase workspace variables, clc->clean command window, close all-> close all currently open figures
clear; clc; close all;
%%
a = 0.09;
K=10;
f = @(t,p) a*p*(1-p./K);
t = 0:2.5:50;
p = -5:0.5:15;

[T,P] = meshgrid(t,p);

dt = t(2) - t(1);
dt2 = dt / 2;
dp = p(2) - p(1);
dp2 = dp / 2;

tmin = t(1) - dt2;
tmax = t(end) + dt2;
pmin = p(1) - dp2;
pmax = p(end) + dp2;

fv = eval(vectorize(f));
yp = feval(fv,T,P);

u = 1./max(1/dt,abs(yp)./dp)*0.35;
v = u .* yp;

quiver(t,p,u,v,0.25,'.k','ShowArrowHead','on');
hold on;
quiver(t,p,-u,-v,0.25,'.k');
hold on;

axis([tmin tmax pmin pmax]);

% Alínea b)
monotonia();

%% Function defs
function monotonia()

    % Variable definition
    a = 0.09; p0 = [-1, 0, 0.1, 5, 10, 15];
    tspan=0:1:50;
    K=10;

    %Tumor growth without drug use 
    for i = 1:6
        [t,p]= ode23s(@(t,p) a*p*(1-p./K),tspan,p0(i),odeset('AbsTol',1e-8,'RelTol',1e-10'));
        plot(t,p,'LineWidth',2)
        hold on
    end

    hold off
    grid; grid minor; xlim([0, 50]);
    title('\textbf{Monotonia da din\^amica de crescimento do tumor}','interpreter','latex','FontSize',32);
    xlabel('\textbf{Tempo, [dias]}','interpreter','latex','FontSize',24); ylabel('\textbf{Crescimento do tumor, [mm$^3$]}','interpreter','latex','FontSize',24);
    legend('','','$V_0 = -1$','$V_0 = 0$', '$V_0 = 0.1$', '$V_0 = 5$','$V_0 = K_T = 10$','$V_0 = 15$','interpreter','latex','FontSize',24, 'NumColumns',3);
end