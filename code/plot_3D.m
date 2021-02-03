%
%  JiZhang.m
%  1D Unsteady
%
%  Created by Yi Qiang Ji Zhang on 19/03/2020.
%  Copyright © 2020 Yi Qiang Ji Zhang. All rights reserved.
%
%  Plot of Temperature Time Distance distribution (using data from JiZhang_plot3D)

clc;
clear;
close all;


%% Plot 3D

% Import data
x = importdata('data/x_data.txt');
y = importdata('data/t_data.txt');
T = importdata('data/temp_data.txt');

% Last instant is not coherent
y(end) = [];
T(end,:) = [];

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
s = surf(X, Y, T);
s.EdgeColor = 'none';
xlim([0 1]);
%ylim([min(y) max(y)+1]);
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Time $\left( \mathrm{days} \right)$");
zlabel("Temperature $\left( \mathrm{K} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
set(gca, 'zticklabel', num2str(get(gca,'ztick')', '%.0f'));
view(75, 20);
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
mesh(X, Y, T);
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
figure(3);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperature - Distance}")
plot(x, T(2191,:), 'b');
xlabel("Distance $\left( \mathrm{m} \right)$");
ylabel("Temperature $\left( ^\circ \mathrm{C} \right)$");
% set(gcf, 'units', 'centimeters', 'position', [18,2,18,18]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
grid on;
grid minor;
box on;
hold off;




