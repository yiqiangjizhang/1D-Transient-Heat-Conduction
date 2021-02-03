%
%  JiZhang.m
%  1D Unsteady
%
%  Created by Yi Qiang Ji Zhang on 19/03/2020.
%  Copyright © 2020 Yi Qiang Ji Zhang. All rights reserved.
%
% Saving data for plot3D_4

clc;
clear;
close all;

%% 1. Data entry

% 1.1 Physical data
rho = 2300.00;                  % [kg/m^3]      % Density
cp = 700.00;                    % [J/kg K]      % Heat at Constant Pressure
k = 1.60;                       % [W/mK]        % Thermal conductivity
h = 8.50;                       % [W/m^2K]      % Heat Transf. coefficient
T0 = 19.00 + 273;               % [K]           % Initial Temperature

% Temperature T1 parameters
A0 = 19.00 + 273;               % [K]
A1 = 5.70;                      % [K]
A2 = 9.10;                      % [K]
w1 = 2*pi/(24*3600);            % [rad/s]
w2 = 2*pi/(365*24*3600);        % [rad/s]

% Domain dimensions
L = 1.5;                         % [m]           % Unitary domain length
H = 1;                          % [m]           % Unitary domain height
W = 1;                          % [m]           % Unitary domain width

% 1.2. Numerical data
N_CV = 250;                     % [units]       % Number of Control Volumes
beta = 0.5;                     % []            % Unsteady Scheme
delta_t = 60;                   % [s]           % Increments of time
t_max = 1*24*60*60;             % [s]           % Maximum time analysis


%% 2. Previous computations

% Number of nodes
N_nodes = N_CV + 2; % [units]

% Distance between nodes
delta_x = L/N_CV; % [m]

% Surface area
surf = H*W; % [m^2]

% Volume control
vol = delta_x*surf; % [m^3]

% Time vector
t = 0:delta_t:t_max; %[s]

% Nodes position vector
nodes_pos = zeros(1,N_nodes); % [m]
for i=1:N_CV % Initial node position and interior nodes position
   nodes_pos(i+1) = (2*i - 1)/2 * delta_x;
end
nodes_pos(end) = L; % Last node position

% Distance vector | dist(i): distance between node 'i'and node 'i+1'
dist = zeros(1,N_CV+1);
for i=1:N_CV+1
    dist(i) = nodes_pos(i+1) - nodes_pos(i);
end


%% 3. Initial Temperatures Map (Physical)

% Temperature map (i = each instant of t) (j = each node)
T = zeros(length(t),N_nodes); % [K]
T(1,:) = T0; % [K]


%% 4. Next instant computation

% External temperature as function of time
T1 = @(t) A0 + A1*sin(w1*t) + A2*sin(w2*t); % [K]

% 4.1. Computing Temperature distribution

for i=1:length(t)-1
    
    % Discretization coefficients definition
    a_W = zeros(1,N_nodes);
    a_E = zeros(1,N_nodes);
    a_P = zeros(1,N_nodes);
    b_P = zeros(1,N_nodes);

    % 4.2. Compute coefficients
    % Initial node (Boundary node)
    a_W(1) = 0;
    a_E(1) = beta*k/(dist(1));
    a_P(1) = beta*h + a_E(1);
    b_P(1) = beta*h*T1(t(i+1)) + (1-beta)*((k/dist(1))*T(i,2) ...
        + h*T1(t(i)) - (h+k/dist(1))*T(i,1));
    
    % Final node (Supposing adiabatic)
    a_W(end) = beta*k/dist(end-1);
    a_E(end) = 0;
    a_P(end) = a_W(end);
    b_P(end) = -(1-beta)*k*((T(end)-T(end-1)))/(dist(end-1));

    % Interior node 
    for j=2:length(nodes_pos)-1
        a_W(j) = beta*k*surf/dist(j-1);
        a_E(j) = beta*k*surf/dist(j);
        a_P(j) = rho*vol*cp/delta_t + a_W(j) + a_E(j);
        b_P(j) = (rho*vol*cp/delta_t - (1-beta)*k*(surf/dist(j-1) + surf/dist(j)))*T(i,j) ...
            + ((1-beta)*k*surf/dist(j-1))*T(i,j-1) + ((1-beta)*k*surf/dist(j))*T(i,j+1);
    end
 
 

    %% 5. Solver using TDMA

    % 5.1. Compute P and Q coefficients
    
    P = zeros(1,N_nodes); % Coefficient P
    Q = zeros(1,N_nodes); % Coefficient Q

    for j=1:N_nodes

        if j==1 % Initial node (Boundary node)
            P(j) = a_E(j)/a_P(j);    
            Q(j) = b_P(j)/a_P(j);

        else % Other nodes
            P(j) = a_E(j)/(a_P(j)-a_W(j)*P(j-1));
            Q(j) = (b_P(j)+a_W(j)*Q(j-1))/(a_P(j)-a_W(j)*P(j-1));  
        end

    end

    % 5.2. Compute temperature

    for j=N_nodes:-1:1 % From the last term to the first
        if j==N_nodes % Last node case
            T(i+1,j)=Q(j);
        else
            T(i+1,j)=P(j)*T(i+1,j+1)+Q(j);      
        end
    end
end

% 6. SAVE DATA TO TEXT FILE

% Just for saving data (Matlab does not allow to do plot in the same script)
% Using N_CV = 1000 is extremely slow (.txt file is over 40 MB) use 100 instead

save data_4 T nodes_pos t;

%% 7. Find the instant of the maximum temperature

maximum = max(max(T1(t)));
[x,y] = find(T1(t) == maximum); % Where y are the indices of the max temperature