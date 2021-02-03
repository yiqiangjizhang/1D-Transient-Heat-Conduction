clear;
close all;
clc;

%% 1. DATA
rho = 2300;
cp = 700;
k = 1.60;
h = 8.50;
T0 = 19 + 273;
A0 = 19 + 273;
A1 = 5.70;
A2 = 9.10;

w1 = 2*pi/(24*3600);
w2 = 2*pi/(365*24*3600);

%% 2. PLOT TEMPERATURE MODES
t1 = 0:60:24*3600;
T_mode0 = A0;
T_mode1 = A1*sin(w1*t1);

t2 = 0:24*3600:365*24*3600;
T_mode2 = A2*sin(w2*t2);

% Mode 1
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature (Mode 1) - Time}");
plot(t1, T_mode1, 'b');
xlabel("Time $\left( \mathrm{hours} \right)$");
ylabel("Temperature $\left( ^\circ\mathrm{C} \right)$");
xlim([0 24*3600]);
xticks([0:2*3600:24*3600]);
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')'/3600, '%.0f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
grid on;
grid minor;
box on;
hold off;

% Mode 2
figure(2);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature (Mode 2) - Temps}");
plot(t2, T_mode2, 'b');
xlabel("Time $\left( \mathrm{days} \right)$");
ylabel("Temperature $\left( ^\circ\mathrm{C} \right)$");
xlim([0 365*24*3600]);
xticks([0:30*24*3600:330*24*3600, 365*24*3600]);
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')'/(24*3600), '%.0f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
grid on;
grid minor;
box on;
hold off;
