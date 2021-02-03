%
%  JiZhang.m
%  1D Unsteady
%
%  Created by Yi Qiang Ji Zhang on 19/03/2020.
%  Copyright © 2020 Yi Qiang Ji Zhang. All rights reserved.
%
%  Plot of Temperature Time Distance distribution (using data from JiZhang_plot3D_6)

clc;
clear;
close all;


%% Plot 3D

% Import data
load('data_6.mat');

x = nodes_pos;
y = t;
z = T;

% Last instant is not coherent
y(end) = [];
z(end,:) = [];

% Change t units to days
y = y/(24*3600);

% For speeding up the graphs (not imprescindible)
% y = y(1:10:end);
% T = T(1:10:end,:);

% Create meshgrid
[X, Y] = meshgrid(x, y);


% Plot 3D of Temperature - Time - Nodes
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature - Time - Distance}")
s = surf(X, Y, z);
s.EdgeColor = 'none';
%xlim([0 1]);
%ylim([min(y) max(y)+1]);
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Time $\left( \mathrm{days} \right)$");
zlabel("Temperature $\left( \mathrm{K} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
% set(gca, 'zticklabel', num2str(get(gca,'ztick')', '%.0f'));
view(50, 13);
grid on;
grid minor;
box on;
hold off;

% Plot of Time - Distance
figure(2);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Time - Distance}")
mesh(X, Y, z);
c = colorbar('TickLabelInterpreter', 'latex');
c.Label.String = sprintf("Temperature $( ^\\circ K )$");
c.Label.Interpreter = 'latex';
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Time $\left( \mathrm{days} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
grid on;
grid minor;
box on;
hold off;

% Plot of Time - Distance
h = figure(3);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature - Distance}")
plot(x, z(2191,:), 'b');
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Temperature $\left( ^\circ \mathrm{C} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
grid on;
grid minor;
box on;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')
hold off;


% Plot 3D of Temperature - Time - Nodes
figure(4);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature - Time - Distance}")
s = surf(X, Y, z);
s.EdgeColor = 'none';
%xlim([0 1]);
%ylim([min(y) max(y)+1]);
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Time $\left( \mathrm{days} \right)$");
zlabel("Temperature $\left( \mathrm{K} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
% set(gca, 'zticklabel', num2str(get(gca,'ztick')', '%.0f'));
view(90, 0);
grid on;
grid minor;
box on;
hold off;

% Plot 3D of Temperature - Time - Nodes
A0 = 19.00 + 273;               % [K]
A1 = 5.70;                      % [K]
A2 = 9.10;                      % [K]
w1 = 2*pi/(24*3600);            % [rad/s]
w2 = 2*pi/(365*24*3600);        % [rad/s]
T_ext = A0 + A1*sin(w1*t) + A2*sin(w2*t); % [K]


figure(5);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature - Time - Distance}")
plot(y,T_ext(1:end-1),'r');
plot(y,z(:,1)','b');
%s.EdgeColor = 'none';
%xlim([0 1]);
%ylim([min(y) max(y)+1]);
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Time $\left( \mathrm{days} \right)$");
%zlabel("Temperature $\left( \mathrm{K} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
% set(gca, 'zticklabel', num2str(get(gca,'ztick')', '%.0f'));
%view(90, 0);
grid on;
grid minor;
box on;
hold off;
