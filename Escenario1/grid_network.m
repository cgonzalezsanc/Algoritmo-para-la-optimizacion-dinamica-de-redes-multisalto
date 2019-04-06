function [links, v_dist] = grid_network(space, n_nodes)
% Creates a grid network

% Check if n_nodes is divisible
max_number = 20;
div = [];
for i=1:max_number
    for j=1:max_number
        if i*j == n_nodes
            div = [div; i j];
        end
    end
end

% We pick the one where the difference between the number of nodes in both
% dimensions is lower
diff_tmp=1e5;
n_nodes_x = 0;
n_nodes_y = 0;
for i=1:size(div,1)
    if abs(div(i,1)-div(i,2)) < diff_tmp
        diff_tmp = abs(div(i,1)-div(i,2));
        n_nodes_x = div(i,1);
        n_nodes_y = div(i,2);
    end
end

% Calculate nodes position in the plane
nodes_pos = [];
dist_x = space(1)/n_nodes_x;
dist_y = space(2)/n_nodes_y;
for i=1:n_nodes_y
    for j=1:n_nodes_x
        nodes_pos = [nodes_pos; j*dist_x, i*dist_y];
    end
end

% Calculate distance matrix
dist_mat = zeros(n_nodes,n_nodes);
for i=1:n_nodes
    for j=1:n_nodes
        dist_mat(i,j) = sqrt((nodes_pos(j,1)-nodes_pos(i,1))^2+(nodes_pos(j,2)-nodes_pos(i,2))^2);
    end
end

% Create adjacency matrix. A node will be only connected to those who are
% above, below, left and right from itself
adj_mat = zeros(n_nodes,n_nodes);
links = [];
for i=1:n_nodes
    for j=1:n_nodes
        if nodes_pos(i,1)-dist_x==nodes_pos(j,1) && nodes_pos(i,2)==nodes_pos(j,2) || ...
           nodes_pos(i,1)+dist_x==nodes_pos(j,1) && nodes_pos(i,2)==nodes_pos(j,2) || ...
           nodes_pos(i,1)==nodes_pos(j,1) && nodes_pos(i,2)-dist_y==nodes_pos(j,2) || ...
           nodes_pos(i,1)==nodes_pos(j,1) && nodes_pos(i,2)+dist_y==nodes_pos(j,2)
            adj_mat(i,j) = 1;
        	links = [links; i, j]; 
        end
    end
end

% % Create adjacency matrix based on nearest nodes
% adj_mat = zeros(n_nodes,n_nodes);
% links = [];
% for i=1:n_nodes
% 	[~, idx] = sort(dist_mat(i,:));
%     neighb = idx(2:linked_nodes+1); % We have to take in consideration the own node
%     adj_mat(i,neighb) = 1;
%     for j=1:linked_nodes
%         links = [links; i neighb(j)];
%     end
% end

% Calculate distance between every link in links
v_dist = zeros(size(links,1),1);
for i=1:length(v_dist)
    v_dist(i) = sqrt((nodes_pos(links(i,1),1)-nodes_pos(links(i,2),1))^2+(nodes_pos(links(i,1),2)-nodes_pos(links(i,2),2))^2);
end

% Plot our net
gplot(adj_mat, nodes_pos, '-*');

end

