%% Programa principal para lanzar las simulaciones
% Escenario 2: Modelo de red distribuido en el que hay varios clusters en
% el plano distribuidos uniformemente con un proceso de Poisson. Dentro de
% cada cluster hay usuarios distribuidos uniformemente cerca de los
% centroides.

clear all; close all;

%% Parámetros variables de la simulación

isDistr = 1;         % Use distributed algorithm
isBatt = 0;          % Use batteries algorithm

% Simulation parameters
showFig = 0;        % Mostrar figuras (asociadas a las baterías)
n_sim = 10;         % Número de simulaciones
storeResults = 0;   % Guardar resultados de la simulación
saveNet = 0;        % Guardar la red creada
loadNet = 0;        % Cargar una red


%% Parámetros fijos de la simulación

space = [100, 100];  % Plane dimensions (meters)
n_prim = 2;          % Número de usuarios primarios
n_interf_links = 2;  % Usuarios que causan interferencia a los primarios
n_sec = 20;          % Number of secondary nodes
n_d2d = 2;           % Number of d2d couples
n_macro = 4;         % Number of macrocells
n_pico = n_macro*4;  % Number of picocells

% Load net
if loadNet
    load('Het n_macro=2 n_pico=8 n_prim=0 n_sec=10 n_d2d=2 space=100 n_interf_links=0.mat'); 
end

% Distributed parameters
eps = 0.1;
d = 0.02;

% Links configuration
if ~exist('links_norm','var')
    [links_norm, links_prim, links_fib, v_dist, v_dist_prim, nodes_pos] = heterogeneous_network(space, n_prim, n_sec, n_d2d, n_macro, n_pico, n_interf_links);
else
    % Add new primary users to the already created network
    if n_prim ~= 0
        [links_prim, v_dist_prim] = add_prim_users_het(space, n_macro, n_pico, n_prim, n_sec, n_d2d, n_interf_links, links_norm, links_fib, links_prim, v_dist_prim, nodes_pos);
    end
end
n_nodes = n_macro+n_pico+n_sec+2*n_d2d;

% Save net
if saveNet
    filename = ['Het n_macro=', num2str(n_macro), ' n_pico=', num2str(n_pico), ' n_prim=', num2str(n_prim), ' n_sec=', num2str(n_sec), ' n_d2d=', num2str(n_d2d), ' space=', num2str(space(1)), ' n_interf_links=', num2str(n_interf_links), '.mat'];
    save(filename,'links_norm','links_fib','links_prim','v_dist','v_dist_prim','nodes_pos');
end

% Flow parameters
dest_flows = [11; 12; 13];  % Destino de cada flujo
src_flows = [14; 15; 16];
n_flows = length(dest_flows);

% Battery and physical parameters
Batt_cap =  10 * ones(n_nodes,1);                % Capacidad máxima de las baterías
Batt_init = Batt_cap/2 + zeros(n_nodes,1);       % Estado inicial de las baterías
P_max =     [10 + zeros(n_macro,1)               % Potencia máxima a transmitir en las macroceldas
             3 + zeros(n_pico,1)                 % En las picoceldas
             1 + zeros(n_sec,1)                  % En los usuarios secundarios
             1 + zeros(2*n_d2d,1)];              % En los D2D
Noises =    1e-10 + zeros(n_nodes,n_nodes);      % Ruido en cada enlace

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
Active_s = zeros(n_nodes,1) + 1;            % Nodos secundarios activos
Average_prim = zeros(n_prim, 1) + 1;        % Posibilidad de transmisión de cada usuario primario
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

sim = struct('n_sec',n_sec,...
             'n_prim',n_prim,...
             'n_macro',n_macro,...
             'n_pico',n_pico,...
             'n_d2d',n_d2d,...
             'n_nodes',n_nodes,...
             'links_norm',links_norm,...
             'links_prim',links_prim,...
             'links_fib',links_fib,...
             'v_dist',v_dist,...
             'v_dist_prim',v_dist_prim,...
             'dest_flows',dest_flows,...
             'src_flows',src_flows,...
             'Batt_cap',Batt_cap,...
             'Batt_init',Batt_init,...
             'P_max',P_max,...
             'Noises',Noises,...
             't',t,...
             'total_duration',total_duration,...
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

mean_rij = zeros(n_nodes,n_nodes,n_flows,n_sim);
a = zeros(n_nodes, n_flows, n_sim);
p = zeros(n_nodes, n_sim);
Opt = zeros(n_sim);
V = zeros(n_sim);
C = zeros(n_sim);
J = zeros(n_sim);

for n = 1:n_sim
    n
    % Ejecución
    [mean_rij(:,:,:,n), a(:,:,n), p(:,n), Opt(n), V(n), J(n), C(n), mean_inter] = sim_alg(sim, isDistr, isBatt, showFig);
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

%% Presentación de resultados

% figure 
% plot(mean_inter(:,1));
% hold on
% plot(Max_interference(:,1), 'g');
% title 'Power recieved by the Primary User'
% xlabel 'Iteration'
% ylabel 'Power recieved (W)'
% legend('Power recieved (W)','Interference limit')
% 
% figure 
% plot(mean_inter(:,2));
% hold on
% plot(Max_interference, 'g');
% title 'Power recieved by the Primary User'
% xlabel 'Iteration'
% ylabel 'Power recieved (W)'
% legend('Power recieved (W)','Interference limit')