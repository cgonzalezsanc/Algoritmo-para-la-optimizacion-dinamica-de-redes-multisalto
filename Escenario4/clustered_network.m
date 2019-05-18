function [links_norm, links_prim, v_dist, v_dist_prim] = clustered_network(space, n_nodes, n_clusters, linked_nodes, n_prim, prim_links)
% Creates a network following a poisson process where users are distributed
% uniformly in space, each one connected to a centroid also distributed
% uniformly.

if n_clusters == 0
    error('n_clusters can not be 0');
end

% Calculate centroid position in the plane
cent_pos = [space(1)*rand(n_clusters, 1), space(2)*rand(n_clusters, 1)];

% Calculate the number of nodes placed in each cluster
n_nodes_clust = repmat(floor(n_nodes/n_clusters), n_clusters,1);
rmd = rem(n_nodes,n_clusters); % Take into account a possible overflow
if rmd ~= 0
    for i=1:rmd
        n = randi(n_clusters); % select randomly the cluster where the node will be added
        n_nodes_clust(n) = n_nodes_clust(n)+1;
    end
end

% Generate the nodes that will be placed in each cluster
links_norm = [];
nodes_pos = [];
max_dist = [space(1), space(2)]./(2*n_clusters);
for i=1:n_clusters
    center = cent_pos(i,:);
    users_pos = [center(1)+(max_dist(1)*(rand(n_nodes_clust(i), 1))-max_dist(1)/2), center(2)+(max_dist(2)*(rand(n_nodes_clust(i), 1))-max_dist(2)/2)];
    nodes_pos = [nodes_pos; users_pos];
%     for j=1:n_nodes_clust(i)
%         nodes_count = nodes_count+1;
%         links_norm = [links_norm; i, nodes_count; nodes_count, i];
%     end
end

% Calculate distance matrix
dist_mat = zeros(n_nodes,n_nodes);
for i=1:n_nodes
    for j=1:n_nodes
        dist_mat(i,j) = sqrt((nodes_pos(j,1)-nodes_pos(i,1))^2+(nodes_pos(j,2)-nodes_pos(i,2))^2);
    end
end

% Create adjacency matrix based on nearest nodes
adj_mat = zeros(n_nodes,n_nodes);
links_norm = [];
for i=1:n_nodes
	[~, idx] = sort(dist_mat(i,:));
    neighb = idx(2:linked_nodes+1); % We have to take in consideration the own node
    adj_mat(i,neighb) = 1;
    for j=1:linked_nodes
        links_norm = [links_norm; i neighb(j)];
    end
end

% Calculate distance between every link in links_norm
v_dist = zeros(size(links_norm,1),1);
for i=1:length(v_dist)
    v_dist(i) = sqrt((nodes_pos(links_norm(i,1),1)-nodes_pos(links_norm(i,2),1))^2+(nodes_pos(links_norm(i,1),2)-nodes_pos(links_norm(i,2),2))^2);
end

% Create adjacency matrix
links_tot = [links_norm];
adj_mat = zeros(n_nodes,n_nodes);
for i=1:length(links_tot)
    adj_mat(links_tot(i,1), links_tot(i,2)) = 1;
end

% Create primary users
prim_pos = [space(1)*rand(n_prim, 1), space(2)*rand(n_prim, 1)];
links_prim = [];
v_dist_prim = [];

if n_prim ~= 0
    % Calculate distance between secondary and primary nodes
    dist_mat_prim = zeros(n_prim,n_nodes);
    for i=1:n_prim
        for j=1:n_nodes
            dist_mat_prim(i,j) = sqrt((nodes_pos(j,1)-prim_pos(i,1))^2+(nodes_pos(j,2)-prim_pos(i,2))^2);
        end
    end

    % Create primary links based on prim_links
    adj_mat_prim = zeros(n_prim,n_nodes);
    for i=1:n_prim
        [dist_sorted, idx] = sort(dist_mat_prim(i,:));
        neighb = idx(1:prim_links);
        adj_mat_prim(i,neighb) = 1;
        for j=1:prim_links
            links_prim = [links_prim; i neighb(j)];
            v_dist_prim = [v_dist_prim; dist_sorted(j)];
        end
    end
end

% Plot our net
figure(1);
%gplot(adj_mat, nodes_pos, '-*');
g = digraph(adj_mat);
plot(g, 'XData', nodes_pos(:,1), 'YData', nodes_pos(:,2));
hold on;
plot(nodes_pos(:,1),nodes_pos(:,2),'r*');
plot(cent_pos(:,1),cent_pos(:,2),'mo');
if n_prim ~= 0
    plot(prim_pos(:,1),prim_pos(:,2),'k^');
    tmp = 0;
    for i=1:n_prim
        for j=1:prim_links
            tmp = tmp+1;
            plot([prim_pos(i,1) nodes_pos(links_prim(tmp,2),1)], [prim_pos(i,2) nodes_pos(links_prim(tmp,2),2)], 'k');
        end
    end
end
axis([0 space(1), 0 space(2)])
legend('Links','Users', 'Centroid', 'Primary users');
grid on;

end

