% Esta función genera la matriz de adyacencia de un conjunto de enlaces

function [Graph] = Create_h_Graph(links_norm, links_fib, nNodes)
Graph = zeros(nNodes, nNodes);
if (max(max(links_norm)) > nNodes)
    Graph = [];
    return;
end;

for i=1:length(links_norm)
    Graph(links_norm(i,1), links_norm(i,2)) = 1;
end;
for i=1:length(links_fib)
    Graph(links_fib(i,1), links_fib(i,2)) = 2;
end;