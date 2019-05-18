function [links_sec, links_prim, v_dist, v_dist_prim, nodes_pos, prim_pos] = poisson_network(space, n_nodes, linked_nodes, n_prim, n_interf_links, nodes_pos, prim_pos)
% Creates a network following a poisson process where users are distributed
% uniformly in space

% Calculate nodes position in the plane
if isempty(nodes_pos)
    nodes_pos = [space(1)*rand(n_nodes, 1), space(2)*rand(n_nodes, 1)];
end

% Calculate distance matrix
dist_mat = zeros(n_nodes,n_nodes);
for i=1:n_nodes
    for j=1:n_nodes
        dist_mat(i,j) = sqrt((nodes_pos(j,1)-nodes_pos(i,1))^2+(nodes_pos(j,2)-nodes_pos(i,2))^2);
    end
end

% Create adjacency matrix and links based on nearest nodes
adj_mat = zeros(n_nodes,n_nodes);
links_sec = [];
for i=1:n_nodes
	[~, idx] = sort(dist_mat(i,:));
    neighb = idx(2:linked_nodes+1); % We have to take in consideration the own node
    adj_mat(i,neighb) = 1;
    for j=1:linked_nodes
        links_sec = [links_sec; i neighb(j)];
    end
end

% Calculate distance between every link in links
v_dist = zeros(size(links_sec,1),1);
for i=1:length(v_dist)
    v_dist(i) = sqrt((nodes_pos(links_sec(i,1),1)-nodes_pos(links_sec(i,2),1))^2+(nodes_pos(links_sec(i,1),2)-nodes_pos(links_sec(i,2),2))^2);
end

% Create primary users
if isempty(prim_pos)
    prim_pos = [space(1)*rand(n_prim, 1), space(2)*rand(n_prim, 1)];
end
links_prim = [];
v_dist_prim = [];

if ~isempty(prim_pos)
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
        neighb = idx(1:n_interf_links);
        adj_mat_prim(i,neighb) = 1;
        for j=1:n_interf_links
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
if n_prim ~= 0
    plot(prim_pos(:,1),prim_pos(:,2),'k^');
    tmp = 0;
    for i=1:n_prim
        for j=1:n_interf_links
            tmp = tmp+1;
            plot([prim_pos(i,1) nodes_pos(links_prim(tmp,2),1)], [prim_pos(i,2) nodes_pos(links_prim(tmp,2),2)], 'k');
        end
    end
end
axis([0 space(1), 0 space(2)])
legend('Links','Users', 'Primary users');
grid on;
saveas(gcf,['Figuras\Poisson n_nodes=', num2str(n_nodes), ' linked_nodes=', num2str(linked_nodes), ' n_prim=', num2str(n_prim), '.png']);

end

