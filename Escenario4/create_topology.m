function [links_norm, links_fib, links_d2d, v_dist] = create_topology(model, n_clusters, space, Nodes, linked_nodes)
% Creates four different topologies:
%   - Poisson distributed network (using uniform distribution)
%   - Clustered distributed network (using uniform distributions)
%   - Grid distributed network
%   - Heterogeneous network with D2D devices

if model == 1
    [links_norm, v_dist] = poisson_network(space, Nodes, linked_nodes);
    links_fib = [];
    links_d2d = [];
elseif model == 2
    [links_norm, links_fib, v_dist] = clustered_network(space, Nodes, n_clusters);
    links_d2d = [];
elseif model == 3
    [links_norm, v_dist] = grid_network(space, Nodes);
    links_fib = [];
    links_d2d = [];
elseif model == 4
    [links_norm, links_fib, links_d2d, v_dist] = heterogeneous_network(space, Nodes, linked_nodes, n_clusters);
end
    
end

