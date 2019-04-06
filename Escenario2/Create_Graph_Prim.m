% Esta funci�n genera la matriz de adyacencia de un conjunto de enlaces

function [Graph, d_Graph] = Create_Graph_Prim(links, n_nodes, n_prim_users, v_dist)

% Inicializaci�n
Graph = zeros(n_prim_users, n_nodes);
d_Graph = zeros(n_prim_users, n_nodes);

% Comprobaci�n de que el n�mero de nodos no sea menor que lo estipulado en
% los enlaces
if (max(max(links)) > n_nodes)
    Graph = [];
    d_Graph = [];
    return;
end

% Se pone un 1 en la posici�n correspondiente a cada enlace y se apunta la
% distancia en la matriz de distancias
for i=1:length(links)
    Graph(links(i,1), links(i,2)) = 1;
    d_Graph(links(i,1), links(i,2)) = v_dist(i);
end

end