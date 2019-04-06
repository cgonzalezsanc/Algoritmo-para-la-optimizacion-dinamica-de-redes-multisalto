%% Programa principal para lanzar las simulaciones
% Escenario 2: Modelo de red distribuido en el que hay varios clusters en
% el plano distribuidos uniformemente con un proceso de Poisson. Dentro de
% cada cluster hay usuarios distribuidos uniformemente cerca de los
% centroides.

clear all; close all;

%% Parámetros variables de la simulación

isDistr = 0;         % Use distributed algorithm
isBatt = 0;          % Use batteries algorithm

% Simulation parameters
showFig = 0;        % Mostrar figuras (asociadas a las baterías)
n_sim = 20;         % Número de simulaciones
storeResults = 0;   % Guardar resultados de la simulación
saveNet = 0;        % Guardar la red creada
loadNet = 0;        % Cargar una red

%% Parámetros fijos de la simulación

n_sec = 10;        % Number of total nodes
linked_nodes = 2;    % Number of directed links per node
n_clusters = 3;      % Number of clusters
space = [100, 100];  % Plane dimensions (meters)
n_prim = 2;          % Número de usuarios primarios
n_interf_links = 2;  % Usuarios que causan interferencia a los primarios

% Load net
if loadNet
    load('Cluster n_clusters=3 n_nodes=10 linked_nodes=2 space=100 n_prim=0 n_interf_links=0.mat'); 
end

% Distributed parameters
eps = 0.1;
d = 0.02;

% Links configuration
if ~exist('links_norm','var')
    [links_norm, links_prim, v_dist, v_dist_prim, nodes_pos] = clustered_network(space, n_sec, n_clusters, linked_nodes, n_prim, n_interf_links);
else
    % Add new primary users to the already created network
    if n_prim ~= 0
        [links_prim, v_dist_prim] = add_prim_users_clust(space, n_sec, n_prim, n_interf_links, links_norm, links_prim, v_dist_prim, nodes_pos);
    end
end
links_fib = [];

% Save net
if saveNet
    filename = ['Cluster n_clusters=', num2str(n_clusters), ' n_nodes=', num2str(n_sec), ' linked_nodes=', num2str(linked_nodes), ' space=', num2str(space(1)), ' n_prim=', num2str(n_prim), ' n_interf_links=', num2str(n_interf_links), '.mat'];
    save(filename,'links_norm','links_prim','v_dist','v_dist_prim','nodes_pos');
end

% Flow parameters
dest_flows = [1; 5; 9];  % Destino de cada flujo
n_flows = length(dest_flows);

% Battery and physical parameters
Batt_cap =  10 * ones(n_sec,1);                % Capacidad máxima de las baterías
Batt_init = Batt_cap/2 + zeros(n_sec,1);       % Estado inicial de las baterías
P_max =     1 + zeros(n_sec,1);                % Potencia máxima a transmitir
Noises =    1e-10 + zeros(n_sec,n_sec);      % Ruido en cada nodo

% Energy distribution parameters
cost_avg = 0.18;                     % Average cost of operator's energy
wind_param = struct('k',2,...        % Shape parameter -> variability of the wind
                    'l',5.65,...     % Scale parameter -> proportional to mean wind speed
                    'p',1.23,...     % Air density
                    'd',1,...        % Rotor blade diameter
                    'Cp', 0.7,...    % Coefficient of performance 
                    'N', 0.7);       % Machinery efficiency 
solar_param = struct('H',1600,...    % Mean solar irradiation (kWh/m^2)
                     'S',1.6,...     % Total solar panel area
                     'r',0.1516,...  % Solar panel efficiency
                     'Pr',0.75);     % Performance ratio

% Time interval
t = 1;                      % Paso de tiempo (en segundos)              
total_duration = 2500;      % Duración de cada simulación

% Primary users parameters
Active_s = zeros(n_sec,1) + 1;            % Nodos secundarios activos
Average_prim = zeros(n_prim, 1) + 1; % Posibilidad de transmisión de cada usuario primario
Max_Interference = 1e-2;                    % Interferencia máxima

% Lagrangian multipliers
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
             'wind_param', wind_param);

%% Lanzamiento de las simulaciones

mean_rij = zeros(n_sec,n_sec,n_flows,n_sim);
a = zeros(n_sec, n_flows, n_sim);
p = zeros(n_sec, n_sim);
Opt = zeros(n_sim);
V = zeros(n_sim);
C = zeros(n_sim);
J = zeros(n_sim);

for n = 1:n_sim
    n
    % Ejecución
    [mean_rij(:,:,:,n), a(:,:,n), p(:,n), Opt(n), V(n), J(n), C(n)] = sim_alg(sim, isDistr, isBatt, showFig);
    x=1;
end

%% Almacenamiento de resultados

if storeResults
    if isBatt
        filename = ['isBatt=', num2str(isBatt), ' Bcap=', num2str(Batt_cap(1)), ' MaxCost=', num2str(max_cost), ' Factor=', num2str(nrg_harv_factor), '.mat'];
        save(filename);
    end
end

% totOpt(tt) = mean(Opt);
% if isDistr ~= 0
%     if storeResults
%         filename = ['isDistr=', num2str(isDistr), ' eps=', num2str(eps), ' d=', num2str(d) ,'.mat'];
%         save(filename, 'Opt', 'totOpt');
%     end    
% end