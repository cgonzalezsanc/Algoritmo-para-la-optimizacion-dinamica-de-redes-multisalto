 %% Programa principal para lanzar las simulaciones
% Escenario 1: Modelo de red distribuido donde la posición de los usuarios
% sigue un proceso de Poisson con distribución uniforme.

clear all; close all;

%% Parámetros que definen el tipo de simulación

% Algoritmos usados
isDistr = 2;         % Use distributed algorithm
isBatt = 0;          % Use batteries algorithm

% Variables para controlar las simulaciones
showFig = 0;         % Mostrar figuras (asociadas a las baterías)
n_sim = 5;          % Número de simulaciones
storeResults = 1;    % Guardar resultados de la simulación
saveNet = 0;         % Guardar la red creada
loadNet = 1;         % Cargar una red
new_links = 0;       % Generar nuevos enlaces

%% Escenario de red

v_n_sec = [10];         % Number of total nodes
v_linked_nodes = [2];    % Number of links of each node
%v_n_prim = [0 0];          % Número de usuarios primarios
% v_Batt_cap = [1 10 1 10 1 10 1 10 1 10 1 10];
% v_cost_avg = [0.1 0.1 0.2 0.2 0.1 0.1 0.2 0.2 0.1 0.1 0.2 0.2];
% v_l = [5.65 5.65 5.65 5.65 11.3 11.3 11.3 11.3 2.32 2.32 2.32 2.32];
% v_H = [1600 1600 1600 1600 2500 2500 2500 2500 800 800 800 800];

for num_scen=1:length(v_n_sec)

% Parámetros del escenario
n_sec = v_n_sec(num_scen);         % Number of total nodes
linked_nodes = v_linked_nodes(num_scen);    % Number of links of each node
%n_prim = v_n_prim(num_scen);          % Número de usuarios primarios
%n_sec = 10;         % Number of total nodes
%linked_nodes = 2;    % Number of links of each node
n_prim = 0;          % Número de usuarios primarios
n_interf_links = 2;   % Usuarios que causan interferencia a los primarios
space = [100, 100];  % Plane dimensions (meters)

% Batt_cap = v_Batt_cap(num_scen);
% cost_avg = v_cost_avg(num_scen);
% l = v_l(num_scen);
% H = v_H(num_scen);
Batt_cap = 10;
cost_avg = 0.18;
l = 5.65;
H = 1600;

% Load net
if loadNet
    %scens = dir('Escenarios');
    %load(strcat('Escenarios\',scens(num_scen+2).name))
    load('02-Poisson n_nodes=10 linked_nodes=2 n_prim=0.mat');
else
    nodes_pos = [];
    prim_pos = [];
end

% Links configuration
if ~exist('links_norm','var')
    [links_norm, links_prim, v_dist, v_dist_prim, nodes_pos, prim_pos] = poisson_network(space, n_sec, linked_nodes, n_prim, n_interf_links, nodes_pos, prim_pos);
else
    if new_links == 1
        [links_norm, links_prim, v_dist, v_dist_prim, nodes_pos, prim_pos] = poisson_network(space, n_sec, linked_nodes, n_prim, n_interf_links, nodes_pos, prim_pos);
    end
    % Add new primary users to the already created network
    if isempty(prim_pos) && n_prim ~= 0
        [prim_pos, links_prim, v_dist_prim] = add_prim_users_poiss(space, n_sec, n_prim, n_interf_links, links_norm, links_prim, v_dist_prim, nodes_pos);
    end
end
links_fib = [];

% Save net
if saveNet
    filename = ['Escenarios\14-Poisson n_nodes=', num2str(n_sec), ' linked_nodes=', num2str(linked_nodes), ' n_prim=', num2str(n_prim), '.mat'];
    save(filename, 'links_norm', 'links_prim', 'v_dist', 'v_dist_prim', 'nodes_pos', 'links_fib', 'n_interf_links', 'space', 'prim_pos');
end

% Time interval
t = 1;                      % Tiempo total(en segundos)              
total_duration = 2500;      % Duración de cada simulación

% Flow parameters
n_flows = 3;
dest_flows = [2; 3; 4];  % Destino de cada flujo

% Primary users parameters
Active_s = zeros(n_sec,1) + 1;            % Nodos secundarios activos
Average_prim = zeros(n_prim, 1) + 1;        % Posibilidad de transmisión de cada usuario primario
Max_Interference = 1e-2;                    % Interferencia máxima

%% Distributed parameters

% CMPA
vEps = [1, 2, 5];%, 1*10^-1, 2*10^-1, 5*10^-1];
vD = [1*10^-2, 2*10^-2, 5*10^-2, 1*10^-1, 2*10^-1, 5*10^-1, 1, 2, 5];
%vEps = 0.001;
%vD = 0.01;

%% Battery and EH parameters

%Batt_cap =  10 * ones(n_sec,1);                % Capacidad máxima de las baterías
Batt_init = Batt_cap/2 + zeros(n_sec,1);       % Estado inicial de las baterías
P_max =     1 + zeros(n_sec,1);                % Potencia máxima a transmitir
Noises =    1e-10 + zeros(n_sec,n_sec);      % Ruido en cada nodo

% Energy distribution parameters
%cost_avg = 0.18;                     % Average cost of operator's energy
wind_param = struct('k',2,...        % Shape parameter -> variability of the wind
                    'l',l,...     % Scale parameter -> proportional to mean wind speed
                    'p',1.23,...     % Air density
                    'd',1,...        % Rotor blade diameter
                    'Cp', 0.7,...    % Coefficient of performance 
                    'N', 0.7);       % Machinery efficiency 
solar_param = struct('H',H,...    % Mean solar irradiation (kWh/m^2)
                     'S',1.6,...     % Total solar panel area
                     'r',0.1516,...  % Solar panel efficiency
                     'Pr',0.75);     % Performance ratio

%% Lagrangian multipliers

if isBatt
    pass_ro = 0.025;                % Paso del multiplicador asociado al enrutado
    pass_pi = 0.25/Batt_cap(1);     % Paso del multiplicador asociado a baterías
    pass_th = 0.05;                 % Paso del multiplicador asociado a interferencia
    ro_init = 0.25;                 % Inicialización de ro
    pi_init = pass_pi*Batt_init;    % Inicialización de pi
    th_init = 0;                    % Inicialización de th
else
    pass_ro = 0.025;                % Paso del multiplicador asociado al enrutado
    pass_pi = 0.025;                % Paso del multiplicador asociado a baterías
    pass_th = 0.05;                 % Paso del multiplicador asociado a interferencia
    ro_init = 0.25;                 % Inicialización de ro
    pi_init = 0.25;                 % Inicialización de pi
    th_init = 0;                    % Inicialización de th
end

%% Generación de la estructura para la simulación

tt = 0;
totOpt = zeros(length(vEps),length(vD));
for ii=1:length(vEps)
for jj=1:length(vD)
    eps = vEps(ii);
    d = vD(jj);
    tt=tt+1;

sim = struct('n_nodes',n_sec,...
             'links_norm',links_norm,...
             'links_prim',links_prim,...
             'links_fib',links_fib,...
             'v_dist',v_dist,...
             'v_dist_prim',v_dist_prim,...
             'dest_flows',dest_flows,...
             'Batt_cap',Batt_cap,...
             'Batt_init',Batt_init,...
             'P_max',P_max,...
             'Noises',Noises,...
             't',t,...
             'total_duration',total_duration,...
             'n_prim',n_prim,...
             'Active_s',Active_s,...
             'Average_prim',Average_prim,...
             'Max_Interference',Max_Interference,...
             'cost_avg',cost_avg,...
             'pass_ro',pass_ro,...
             'pass_pi',pass_pi,...
             'pass_th',pass_th,...
             'ro_init',ro_init,...
             'pi_init',pi_init,...
             'th_init',th_init,...
             'eps',eps,...
             'd',d,...
             'solar_param',solar_param,...
             'wind_param', wind_param,...
             'num_scen', num_scen,...
             'linked_nodes',linked_nodes);

%% Lanzamiento de las simulaciones

rij = zeros(n_sec,n_sec,n_flows,n_sim);
a = zeros(n_sec, n_flows, n_sim);
p = zeros(n_sec, n_sim);
Opt = zeros(n_sim,1);
V = zeros(n_sim,1);
C = zeros(n_sim,1);
J = zeros(n_sim,1);
e_out = zeros(n_sec,n_sim);
e_in = zeros(n_sec,n_sim);
e_op = zeros(n_sec,n_sim);
batt = zeros(n_sec,n_sim);

for n = 1:n_sim
    n
    % Ejecución
    [rij(:,:,:,n), a(:,:,n), p(:,n), Opt(n), V(n), J(n), C(n), e_out(:,n), e_in(:,n), e_op(:,n), batt(:,n)] = sim_alg(sim, isDistr, isBatt, showFig, n);
end

mean_rij = mean(rij,4);
mean_a = mean(a,3);
mean_p = mean(p,2);
mean_Opt = mean(Opt);
mean_V = mean(V);
mean_J = mean(J);
mean_C = mean(C);
mean_e_out = mean(e_out,2);
mean_e_in = mean(e_in,2);
mean_e_op = mean(e_op,2);
mean_batt = mean(batt,2);

%% Almacenamiento de resultados

totOpt(ii,jj) = mean(Opt);
if isDistr ~= 0
    if storeResults
        filename = ['Resultados\',num2str(tt+108),'-n_sec=', num2str(n_sec), ' eps=', num2str(eps), ' d=', num2str(d) ,'.mat'];
        save(filename, 'Opt', 'totOpt');
    end    
end

% totOpt(tt) = mean(Opt);
% if isDistr ~= 0
%     if storeResults
%         filename = ['isDistr=', num2str(isDistr), ' eps=', num2str(eps), ' d=', num2str(d) ,'.mat'];
%         save(filename, 'Opt', 'totOpt');
%     end    
% end
end
end
end