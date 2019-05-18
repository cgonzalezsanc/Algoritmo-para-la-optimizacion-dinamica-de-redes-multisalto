function [links_prim, v_dist_prim] = add_prim_users_het(space, n_micro, n_pico, n_prim, n_sec, n_d2d, n_interf_links, links_norm, links_fib, links_prim, v_dist_prim, nodes_pos)
% Function that adds new primary user to an existing network

% Create primary users
prim_pos = [space(1)*rand(n_prim, 1), space(2)*rand(n_prim, 1)];

if n_prim ~= 0
    % Calculate distance between secondary and primary nodes
    dist_mat_prim = zeros(n_prim,n_d2d+n_sec);
    for i=1:n_prim
        for j=n_pico+n_micro+1:n_pico+n_micro+n_d2d+n_sec
            dist_mat_prim(i,j-n_pico-n_micro) = sqrt((nodes_pos(j,1)-prim_pos(i,1))^2+(nodes_pos(j,2)-prim_pos(i,2))^2);
        end
    end

    % Create primary links based on prim_links
    adj_mat_prim = zeros(n_prim,n_d2d+n_sec);
    for i=1:n_prim
        [dist_sorted, idx] = sort(dist_mat_prim(i,:));
        neighb = idx(1:n_interf_links);
        adj_mat_prim(i,neighb) = 1;
        for j=1:n_interf_links
            links_prim = [links_prim; i neighb(j)+n_pico+n_micro];
            v_dist_prim = [v_dist_prim; dist_sorted(j)];
        end
    end
end

% Split the position of the different type of nodes
micros_pos = nodes_pos(1:n_micro,:);
picos_pos = nodes_pos(n_micro+1:n_micro+n_pico,:);
sec_users_pos = nodes_pos(n_micro+n_pico+1:n_micro+n_pico+n_sec,:);
d2d_users_pos = nodes_pos(n_micro+n_pico+n_sec+1:end,:);

% Create adjacency matrix for fiber links
n_bs = n_micro+n_pico;
adj_mat_fib = zeros(n_bs,n_bs);
for i=1:length(links_fib)
    adj_mat_fib(links_fib(i,1), links_fib(i,2)) = 1;
end
bs_pos = [micros_pos;picos_pos];

% Create adjacency matrix for secondary links
n_wless = n_sec+2*n_d2d+n_bs;
adj_mat_wless = zeros(n_wless,n_wless);
for i=1:length(links_norm)
    adj_mat_wless(links_norm(i,1), links_norm(i,2)) = 1;
end
wless_pos = [bs_pos;sec_users_pos;d2d_users_pos];

% Plot our net
figure(1);
%gplot(adj_mat_wless, wless_pos, 'b-');
g_wless = digraph(adj_mat_wless);
plot(g_wless, 'XData', wless_pos(:,1), 'YData', wless_pos(:,2));
hold on;
gplot(adj_mat_fib, bs_pos, 'r-');
%g_fib = digraph(adj_mat_fib);
%plot(g_fib, 'XData', bs_pos(:,1), 'YData', bs_pos(:,2));
plot(micros_pos(:,1),micros_pos(:,2),'ro');
plot(picos_pos(:,1),picos_pos(:,2),'gp');
plot(sec_users_pos(:,1),sec_users_pos(:,2),'k*');
plot(d2d_users_pos(:,1),d2d_users_pos(:,2),'m^');
if n_prim ~= 0
    plot(prim_pos(:,1),prim_pos(:,2),'rx');
    tmp = 0;
    for i=1:n_prim
        for j=1:n_interf_links
            tmp = tmp+1;
            plot([prim_pos(i,1) nodes_pos(links_prim(tmp,2),1)], [prim_pos(i,2) nodes_pos(links_prim(tmp,2),2)], 'r:');
        end
    end
end
axis([0 space(1), 0 space(2)])
legend('Wireless links','Fiber links','Microcells', 'Picocells', 'Secondary users', 'D2D users', 'Primary users', 'Interference links');
grid on;

end

