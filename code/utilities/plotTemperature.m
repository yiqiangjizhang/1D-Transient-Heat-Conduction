clear;
close all;
clc;

%% 1. DATA
T_scale = 1.001;
T0 = 0;
L = 1;
N = 1e2;

x = linspace(-L, L, N);
y = linspace(-L, L, N);

[X,Y] = meshgrid(x, y);
surface = T_scale*cos((2*pi/L)*X).*sin((2*pi/L)*Y) + T0;

%% 2. TEMPERATURE COMPUTATION
T = zeros(N, N);
for j = 1:N
    T(j,:) = T_scale*cos((2*pi/L)*x)*sin((2*pi/L)*y(j)) + T0;
end

%% 3. PLOT
null = zeros(1,N);

figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Distribuci\'o de Temperatures}");
% Temperature plot
    s = surf(X, Y, surface);
    s.EdgeColor = 'none';
% Right colorbar
    c = colorbar('TickLabelInterpreter', 'latex');
    c.Label.String = sprintf("Temperatura $( ^\\circ C )$");
    c.Label.Interpreter = 'latex';
    c.Ticks = [-1.2:0.2:1.2]; %#ok<NBRAK>
    realTicks = get(c,'xtick');
    colorbarTicks = cell(1,length(realTicks));
    for i = 1:length(get(c,'xtick'))
        colorbarTicks(i) = {sprintf("%.2f", realTicks(i))};
    end
    c.TickLabels = colorbarTicks;
    c.Limits = [-1 1];
% Plot labels
xlabel("$x \ \left( \mathrm{m} \right)$");
ylabel("$y \ \left( \mathrm{m} \right)$");
zlabel("$T \ \left( ^\circ \mathrm{C} \right)$");
zlim([-1 1]);
set(gcf, 'units', 'centimeters', 'position', [18,1,20,18]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
view(45, 30);
grid on;
grid minor;
box on;
hold off;