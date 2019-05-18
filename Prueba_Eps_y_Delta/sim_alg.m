function [mean_rij, a, p_mean, Opt, V, J, C, mean_e_out, mean_e_in, mean_e_op, mean_rec_batt] = sim_alg (sim, isDistr, isBatt, showFig, num_sim)
% Parametros de entrada


%%

%====================================
%===================================
% PHASE I
%
% Definition of the Network
%====================================
%===================================

%------------------------------------
% Links and Graph definition
%------------------------------------

% Links del medio móvil
links_norm = sim.links_norm;
% Links conectados por medio cableado
links_fib = sim.links_fib;
% Links con los usuarios primarios
links_prim = sim.links_prim;
% Vector de distancias
v_dist = sim.v_dist;
% Vector de distancias con usuarios primarios
v_dist_prim = sim.v_dist_prim;

% Conjunto de todos los links
links = [links_norm; links_fib];

% Número de nodos
n_nodes = sim.n_nodes;
% Creación de la matriz de adyacencia del grafo y de la matriz de
% distancias
[Graph, d_Graph] = Create_Graph(links, n_nodes, v_dist);

% El número del flujo indica su destino
dest_flows = sim.dest_flows;

%------------------------------------
% Physical parameters
%------------------------------------

% Noise parameters
Sigma=200;
Error=0;
Error_G=0;
window=0;

% Maxima interferencia a los PU
Max_Interference = sim.Max_Interference;
% Capacidad de la batería de cada nodo
Batt_cap = sim.Batt_cap;
% Potencia máxima a transmitir por enlace y flujo
P_max = sim.P_max;


%------------------------------------
% Definition of the Flows
%------------------------------------

N_Flows = length(dest_flows);

%------------------------------------
% Special characterictics of some Nodes
%------------------------------------

Noises = sim.Noises;
Noises = Noises.*Graph;

min_tasa = zeros(n_nodes, N_Flows);

%------------------------------------
% Definition of the primary Users
%------------------------------------
n_prim_users = sim.n_prim;
Active_s = sim.Active_s;
Average_prim = sim.Average_prim;

%------------------------------------
% Cost of operator energy
%------------------------------------

cost_avg = sim.cost_avg;

%------------------------------------
% Set energy harvested by user 
%------------------------------------

solar_param = sim.solar_param;
wind_param = sim.wind_param;

%------------------------------------
% Distributed parameters
%------------------------------------

eps = sim.eps;
d = sim.d;

%====================================
%===================================
% PHASE II
%
% Ergodic parameters initialization
%===================================
%===================================

%------------------------------------
% Lagrange multipliers
%------------------------------------
pass_ro = sim.pass_ro;
pass_pi = sim.pass_pi;
pass_th = sim.pass_th;

ro_init = sim.ro_init;
pi_init = sim.pi_init;
th_init = sim.th_init;

%------------------------------------
% Utility functions are

% V(a) = ln(0.1 + a)
% J(p) =  p
% C(p_op) = nrg_cost*rec_pi_op
%------------------------------------
%%
%====================================
%===================================
% PHASE III
%
% Initialization of the algorithm
%===================================
%===================================

%------------------------------------
% Channels definition
%------------------------------------
t = sim.t;
total_duration = sim.total_duration;
it = t/total_duration;

hGraph = Create_h_Graph(links_norm, links_fib, n_nodes);

h = generate_channels(hGraph, t, it, Sigma, d_Graph);
h_error = generate_channels(Graph, t, it, Error, d_Graph);

h_est = h + h_error;

%------------------------------------
% Channel SU-PU definition
%------------------------------------

% Creación de la matriz de adyacencia del grafo y de la matriz de
% distancias con los usuarios primarios
[Graph_Prim, d_Graph_Prim] = Create_Graph_Prim(links_prim, n_nodes, n_prim_users, v_dist_prim);

G = generate_p_channels(Graph_Prim, t, it, Sigma, n_nodes, n_prim_users, d_Graph_Prim);
G_Error = generate_p_channels(Graph_Prim, t, it, Error_G, n_nodes, n_prim_users, d_Graph_Prim);

G_est = G + G_Error;

%------------------------------------
% Formation of the conflict graph 
%------------------------------------

% Si no es distribuido se genera directamente el conjunto de MIS.
% Si es distribuido se crea el grafo de conflicto.
if ~isDistr
    [Sets, N_Sets] = Independent_Sets(links_norm, links_fib, Graph);
else
    conflict_graph = Create_conflict_graph(links_norm, links_fib, Graph);
end

%------------------------------------
% Energy calculations
%------------------------------------

nrg_cost = get_energy_cost(total_duration, cost_avg);
nrg_harv = get_energy_harvested(n_nodes, total_duration, solar_param, wind_param, it);

%------------------------------------
% Ergodic parameters and trafic data
%------------------------------------
ro = zeros (n_nodes, N_Flows) + ro_init;

for f = 1:N_Flows
    ro(dest_flows(f), f) = 0;
end;
pi = zeros (n_nodes, 1) + pi_init;
pi_e = zeros (n_nodes, 1);
th = zeros (n_prim_users, 1) + th_init/n_prim_users;

a = zeros (n_nodes, N_Flows);
p = zeros (n_nodes,1);

%------------------------------------
% Routing parameters
%------------------------------------
mean_rij = zeros (length(Graph), length(Graph), N_Flows);
record_rij_nodes_flows = zeros(total_duration, n_nodes, n_nodes, N_Flows);
rij_nodes_flows = zeros(total_duration, n_nodes, n_nodes, N_Flows);

%------------------------------------
% Interference parameters
%------------------------------------
acum_inter = zeros (total_duration, n_prim_users);
mean_inter = zeros (total_duration, n_prim_users);

%------------------------------------
% Physical parameters
%------------------------------------
rec_batt = zeros(n_nodes,total_duration+1); % state of the battery
rec_batt(:,1) = sim.Batt_init;              % battery initialization
rec_pi = zeros(n_nodes, total_duration);    % record of pi multiplier
rec_pi(:,1) = pi;                           % initialization of pi multiplier
rec_p = zeros(n_nodes, total_duration);     % power transmitted by user i
rec_e_op = zeros(n_nodes, total_duration);  % power bought by user i
rec_e_out = zeros(n_nodes, total_duration); % power taken out from the battery

%====================================
%%
%===================================
%===================================
% Begining of the algorithm
%===================================
%===================================


phi1=0;

for i = 1:total_duration
   
    if mod(i, 500)==0
        i
    end
    
    %-----------------
    %Init
    %-----------------
    
    %Physycal layer values
    rij = zeros (n_nodes, n_nodes, N_Flows);
    cij_h = zeros (n_nodes, n_nodes);
    pij_h = zeros (n_nodes, n_nodes);
    phi = zeros (n_nodes, n_nodes);
    pot_tot = zeros(n_nodes, 1);
    
    % Harvesting values
    e_op = zeros(n_nodes,1);
    e_in = zeros(n_nodes,1);
    e_out = zeros(n_nodes,1);

    % Set values
    wij = zeros(n_nodes, n_nodes);
    
    %Ergodic definitions
    ro_ij = zeros (n_nodes, n_nodes);
    sf = zeros(n_nodes, n_nodes);       
    pi_tot = pi + pi_e;
    
    %-----------------
    %First part. Getting phi
    %-----------------
    
    for x1=1:n_nodes
        for x2=1:n_nodes
            if Graph(x1,x2) == 1
                
                %Get ro_ij and best flow for each link
                diff_ro = zeros(N_Flows, 1);
                for f = 1:N_Flows
                    diff_ro(f) = ro(x1,f) - ro(x2,f);
                    if (dest_flows(f) == x1)
                        diff_ro(f) = -999;
                    end;
                end; 
                [bst, loc_bst] = max(diff_ro);
                ro_ij(x1,x2) = Graph(x1,x2).*bst;
                sf(x1,x2) = Graph(x1,x2).*loc_bst;
                
                %Get phi, pij_h and cij_h
                [phi(x1,x2), pij_h(x1,x2), cij_h(x1,x2)] = Calcular_param(h_est(i,x1,x2), ro_ij(x1,x2), pi_tot(x1,1), P_max(x1), Noises(x1,x2));
            end;
        end;
    end;
    
    %-----------------
    %Select the best Set to transmit
    %-----------------
   
    if isDistr == 1
        [Sets, N_Sets] = Ind_Set_MWDGA(conflict_graph, phi, links_norm, links_fib);
    elseif isDistr == 2
        [Sets, N_Sets] = Ind_Set_CMPA(conflict_graph, phi, links_norm, links_fib, eps, d);
        %[Sets1, N_Sets1] = Ind_Set_MaxProd(conflict_graph, phi, links_norm, links_fib);
        %[Sets2, N_Sets2] = Ind_Set_MWDGA(conflict_graph, phi, links_norm, links_fib);
    end
    
    phi_s = zeros(N_Sets, 1);
    w_s = zeros(N_Sets, 1);
    if ~isDistr
        %Formation of the phi parameter for each Independet Set
        for is=1:N_Sets
            lk = 1;
            while (Sets(is,lk) > 0)
                x1 = links(Sets(is,lk), 1);
                x2 = links(Sets(is,lk), 2);

                phi_s(is, 1) = phi_s(is, 1) + phi(x1,x2);
                lk = lk + 1;
            end;
        end;
    else
       for lk = 1:length(Sets)
           x1 = links(Sets(lk), 1);
           x2 = links(Sets(lk), 2);
           phi_s = phi_s+phi(x1,x2);
       end
    end
    phi1 = phi1+phi_s;
    
    %Select the best Independent set
    [phisel, selected] = min(phi_s);
    
    if length(find(phi_s == min(phi_s))) > 1
        %Draw
        draws = length(find(phi_s == min(phi_s)));
        to_select = find(phi_s == min(phi_s));
        random_number = ceil ( draws * rand(1));
        selected = to_select(random_number);
    end;
    
    w_s(selected) = 1;
    if ~isDistr
        for is=1:N_Sets
            lk = 1;
            while (Sets(is,lk) > 0)
                x1 = links(Sets(is,lk), 1);
                x2 = links(Sets(is,lk), 2);
                wij(x1,x2) = wij(x1,x2) + w_s(is,1);  
                lk = lk + 1;
            end;
        end;
    else
       for lk = 1:length(Sets)
           x1 = links(Sets(lk), 1);
           x2 = links(Sets(lk), 2);
           wij(x1,x2) = wij(x1,x2) + w_s;
       end        
    end
    
    %-----------------
    %Transmission
    %-----------------
    
    %Definition of the power and routing values for the transmission
    pij = pij_h .* wij;
    cij = cij_h .* wij;
    
    for m = 1:n_nodes
        for n=1:n_nodes
            if (wij(m,n) == 1)
                bst = sf(m,n);
                rij(m,n,bst) = cij(m,n);
            end;
        end;
    end;
 
    %Total power transmitted by each node
    for m = 1:n_nodes
        pot_tot(m,1) = sum(pij(m,:),2);
        rec_p(m,i) = pot_tot(m,1);
    end;
    
    %-----------------
    % Harvesting algorithm
    %-----------------
    
    if isBatt
        for m = 1:n_nodes
           % Collect energy
           [e_op(m), e_out(m), e_in(m)] = harvesting_algorithm(pot_tot(m), nrg_harv(m,i), nrg_cost(i), Batt_cap, rec_batt(m,i), pi(m));
           % Update battery status
           rec_batt(m,i+1) = rec_batt(m,i) + e_in(m) - e_out(m);
           rec_e_op(m,i) = e_op(m);
           rec_e_out(m,i) = e_out(m);
        end
    end
    
    %-----------------
    %Update of Ergodic parameters
    %-----------------
    
    %Update the average traffic routed
    rij_nodes_flows(i,:,:,:) = rij;
    for m = 1:length(Graph)
        for n=1:length(Graph)
            for f = 1:N_Flows
                mean_rij(m,n,f) = mean(rij_nodes_flows(:,m,n,f));
                record_rij_nodes_flows(i,m,n,f) = mean_rij(m,n,f);
            end;
        end;
    end;

    %Get interference caused to every PU update interference parametrers
    for k = 1:n_prim_users
        A_Interference = 0;
        A_Interference_r = 0;
         for m = 1:n_nodes
             if abs(G_est(i,k,m))^2 > 0
                  A_Interference = A_Interference + Average_prim(k,1) * Active_s(m) * pot_tot(m) * abs(G_est(i,k,m))^2;
                  A_Interference_r = A_Interference_r + Average_prim(k,1) * Active_s(m) * pot_tot(m) * abs(G(i,k,m))^2;
             end;
         end;
         acum_inter (i,k) = A_Interference_r;
         th(k, 1) = max(th(k,1) + pass_th * (A_Interference  - Max_Interference),0);
    end;
    
    %Update of ergodic parameters
    for m = 1:n_nodes
        
        %parameter ro
        for f = 1:N_Flows
            a(m,f) = max(1 /(0.05 +  ro(m,f)) - 1/10.05,0);
            
            a(dest_flows(f), f) = 0;
                  
            if a(m,f) < min_tasa(m,f)
                a(m,f) = min_tasa(m,f);
            end;
            
            ro(m,f) = min(ro(m,f) + pass_ro*( a(m,f) + sum(mean_rij(:,m,f)) - sum(mean_rij(m,:,f))) , 10);
            if ro(m,f) < 0
                ro(m,f) = 0;
            end;
        end;
        
        %parameter pi
        p(m,1) = pi(m,1);
        if ~isBatt
            pi(m,1) = max(pi(m,1) + pass_pi * (sum(pij(m,:)) - p(m,1)),0)+1e-6*rand(1,1);
        else
            C_max = pass_pi*Batt_cap;
            pi(m,1) = max(C_max-pass_pi*rec_batt(m,i+1), 0)+1e-6*rand(1,1);
        end
        rec_pi(m,i+1) = pi(m,1);
        
        %parameter pi_e
        pi_e(m,1) = 0;
        for k = 1:n_prim_users
            if abs(G_est(i,k,m))^2 > 0
                pi_e(m,1) = pi_e(m,1) + mean(abs(G(:,k,m)).^2) * (th(k,1));
            end;
        end;
    end;
    
    for f = 1:N_Flows
        ro(dest_flows(f), f) = 0;
    end;
    
    %Record of the average interference recieved   
    for k = 1:n_prim_users
          if window == 0
            mean_inter(i,k) = mean(acum_inter(10:i,k));
          else 
              if i > 25
                  mean_inter(i,k) = mean(acum_inter(i-25:i,k));
              else 
                  mean_inter(i,k) = mean(acum_inter(10:i,k));
              end;
          end;
    end;   
  
           
end;
%%

p_mean = mean(rec_p.');

if showFig && num_sim == 15
    for i=1:n_nodes
        close all;

        t=1:2000;

        figure(2);
        plot(t, nrg_cost(t), 'b'); 
        hold on; 
        plot(t, rec_pi(i,t), 'r');
        titulo = ['Comparativa entre el coste de comprar potencia y pi en el nodo ', num2str(i)];
        title(titulo);
        legend('Coste','\Pi');
        saveas(gcf,['Figuras\Baterias\Batt_cap=',num2str(Batt_cap),' cost_avg=', num2str(cost_avg), ' l=', num2str(wind_param.l), ' H=', num2str(solar_param.H), '\' ,titulo,'.png']);

        figure(3);
        plot(t, rec_e_op(i,t)-rec_e_out(i,t), 'b'); 
        titulo = ['Potencia comprada - Potencia sacada de la batería en el nodo ', num2str(i)];
        title(titulo);
        ylim([-1 1]);
        saveas(gcf,['Figuras\Baterias\Batt_cap=',num2str(Batt_cap),' cost_avg=', num2str(cost_avg), ' l=', num2str(wind_param.l), ' H=', num2str(solar_param.H), '\' ,titulo,'.png']);
        
        figure(4);
        plot(t,rec_batt(i,t), 'b');
        titulo = ['Estado de la batería en el nodo ', num2str(i)];
        title(titulo);
        saveas(gcf,['Figuras\Baterias\Batt_cap=',num2str(Batt_cap),' cost_avg=', num2str(cost_avg), ' l=', num2str(wind_param.l), ' H=', num2str(solar_param.H), '\' ,titulo,'.png']);
    end
end

V = 0;
J = 0;
C = 0;
for m =1:n_nodes
    for f = 1:N_Flows
        if (a(m,f) >0 )
        V = V + log(0.1 + a(m,f));
        end;
    end;
    J = J + mean(rec_p(m,:)).^2;
    C = C + mean(nrg_cost.*rec_e_op(m,:));
end;

mean_e_out = mean(e_out,2);
mean_e_in = mean(e_in,2);
mean_e_op = mean(e_op,2);
mean_rec_batt = mean(rec_batt,2);

Opt = V - J - C

% figure 
% plot(mean_inter(:,1));
% hold on
% plot(int(:,1), 'g');
% title 'Power recieved by the Primary User'
% xlabel 'Iteration'
% ylabel 'Power recieved (W)'
% legend('Power recieved (W)','Interference limit')
% 
% figure 
% plot(mean_inter(:,2));
% hold on
% plot(int(:,1), 'g');
% title 'Power recieved by the Primary User'
% xlabel 'Iteration'
% ylabel 'Power recieved (W)'
% legend('Power recieved (W)','Interference limit')
%  


record_rij_flows = zeros(total_duration, N_Flows);

for f=1:N_Flows
    for i=1:total_duration
       record_rij_flows(i,f) = sum(sum(record_rij_nodes_flows(i,:,:,f)));
    end;
end;

sf = sum(record_rij_flows(total_duration,:));
interfer = mean_inter(total_duration,:);

%filename = ['Resultados\', num2str(sim.num_scen+22), '-Poisson n_nodes=',num2str(n_nodes),' linked_nodes=', num2str(sim.linked_nodes), '-Total.mat'];
%save(filename);

end