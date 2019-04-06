function [links_prim, v_dist_prim] = add_prim_users_poiss(space, n_nodes, n_prim, n_interf_links, links_norm, links_prim, v_dist_prim, nodes_pos)
% Function that adds new primary user to an existing network

% Create primary users
prim_pos = [space(1)*rand(n_prim, 1), space(2)*rand(n_prim, 1)];

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
        neighb = idx(1:n_interf_links);
        adj_mat_prim(i,neighb) = 1;
        for j=1:n_interf_links
            links_prim = [links_prim; i neighb(j)];
            v_dist_prim = [v_dist_prim; dist_sorted(j)];
        end
    end
end

% Create adjacency matrix based on links_norm
adj_mat = zeros(n_nodes,n_nodes);
for i=1:length(links_norm)
    adj_mat(links_norm(i,1), links_norm(i,2)) = 1;
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

end

