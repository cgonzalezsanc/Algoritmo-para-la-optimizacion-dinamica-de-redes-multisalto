function [links_sec, links_prim, v_dist, v_dist_prim, nodes_pos] = grid_network(space, n_nodes, n_prim, prim_links)
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

% We pick the division where the difference between the number of nodes in both
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

% Create adjacency matrix. A node will be only connected to those who are
% above, below, left and right from itself
adj_mat = zeros(n_nodes,n_nodes);
links_sec = [];
for i=1:n_nodes
    for j=1:n_nodes
        if check_conditions(nodes_pos(i,:),nodes_pos(j,:),dist_x,dist_y)
            adj_mat(i,j) = 1;
        	links_sec = [links_sec; i, j]; 
        end
    end
end

% Calculate distance between every link in links
v_dist = zeros(size(links_sec,1),1);
for i=1:length(v_dist)
    v_dist(i) = sqrt((nodes_pos(links_sec(i,1),1)-nodes_pos(links_sec(i,2),1))^2+(nodes_pos(links_sec(i,1),2)-nodes_pos(links_sec(i,2),2))^2);
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
axis([0 space(1)+dist_x, 0 space(2)+dist_y])
legend('Links','Users', 'Primary users');
grid on;

end

function bool = check_conditions(node_1,node_2,dist_x,dist_y)
% We introduce a threshold so we are not vulnerable to small decimals
bool = (abs(node_1(1)-dist_x-node_2(1)) < 1e-5 && ...
        node_1(2)==node_2(2)) || ...
       (abs(node_1(1)+dist_x-node_2(1)) < 1e-5  && ...
        node_1(2)==node_2(2)) || ...
       (node_1(1)==node_2(1) && ...
        abs(node_1(2)-dist_y-node_2(2)) < 1e-5 ) || ...
       (node_1(1)==node_2(1) && ...
        abs(node_1(2)+dist_y-node_2(2)) < 1e-5 );
end
