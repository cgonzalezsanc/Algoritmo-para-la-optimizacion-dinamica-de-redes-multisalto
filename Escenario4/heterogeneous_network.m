function [links_norm, links_prim, links_fib, v_dist, v_dist_prim, nodes_pos] = heterogeneous_network(space, n_prim, n_sec, n_d2d, n_macro, n_pico, prim_links)
% Creates an heterogeneous network which includes:
%   - macrocells (fixed distribution)
%   - Picocells (fixed distribution)
%   - Secondary users (uniformly distributed)
%   - D2D users (uniformly distributed)
%   - Primary users (uniformly distributed)

% -----------------------------------
% Calculate macrocell position
if n_macro == 0
    error('n_macro can not be 0');
end

% Check if n_macro is divisible
max_number = 20;
div = [];
for i=1:max_number
    for j=1:max_number
        if i*j == n_macro
            div = [div; i j];
        end
    end
end

% We pick the division where the difference between the number of macros in both
% dimensions is lower
diff_tmp = 1e5;
n_macro_x = 0;
n_macro_y = 0;
for i=1:size(div,1)
    if abs(div(i,1)-div(i,2)) < diff_tmp
        diff_tmp = abs(div(i,1)-div(i,2));
        n_macro_x = div(i,1);
        n_macro_y = div(i,2);
    end
end

% Calculate macrocell position in the plane
macros_pos = [];
dist_macro_x = space(1)/n_macro_x;
start_macro_x = dist_macro_x/2;
dist_macro_y = space(2)/n_macro_y;
start_macro_y = dist_macro_y/2;
for i=1:n_macro_y
    for j=1:n_macro_x
        macros_pos = [macros_pos; start_macro_x+(j-1)*dist_macro_x, start_macro_y+(i-1)*dist_macro_y];
    end
end
nodes_pos = macros_pos;

% macrocells are connected by fiber
links_fib = [];
for i=1:n_macro
    for j=1:n_macro
        if i ~= j
            links_fib = [links_fib; i j];
        end
    end
end

% Calculate picocell position. For simplicity, we assume there are 4 in
% every macrocell.
dist_pico_x = dist_macro_x/2;
dist_pico_y = dist_macro_y/2;
picos_pos = [];
for i=1:n_macro
    picos_pos = [picos_pos
                 macros_pos(i,1)-dist_pico_x/2, macros_pos(i,2)-dist_pico_y/2
                 macros_pos(i,1)+dist_pico_x/2, macros_pos(i,2)-dist_pico_y/2
                 macros_pos(i,1)-dist_pico_x/2, macros_pos(i,2)+dist_pico_y/2
                 macros_pos(i,1)+dist_pico_x/2, macros_pos(i,2)+dist_pico_y/2];
end
nodes_pos = [nodes_pos; picos_pos];

% Picocells are connected to macrocells by fiber
for i=1:n_pico
    idx_pico = n_macro+i;
    idx_macro = floor((i-1)/(n_pico/n_macro))+1;
    links_fib = [links_fib; idx_macro, idx_pico; idx_pico, idx_macro];
end

% Generate secondary users uniformly distributed over the plane
sec_users_pos = [space(1)*rand(n_sec, 1), space(2)*rand(n_sec, 1)];
nodes_pos = [nodes_pos; sec_users_pos];

% Every secondary user has to be connected to the nearest picocell and the
% nearest macrocell
links_norm = [];
for i=1:n_sec
    n_user = i+n_pico+n_macro;
    % Get distance to macrocells
    dist_to_macro = sqrt((sec_users_pos(i,1)-macros_pos(:,1)).^2+(sec_users_pos(i,2)-macros_pos(:,2)).^2);
    % Get nearest macro
    [~, idx] = sort(dist_to_macro);
    nearest_macro = idx(1);
    % Get distance to picocells
    dist_to_pico = sqrt((sec_users_pos(i,1)-picos_pos(:,1)).^2+(sec_users_pos(i,2)-picos_pos(:,2)).^2);
    % Get nearest pico
    [~, idx] = sort(dist_to_pico);
    nearest_pico = idx(1)+n_macro;
    % Link the user to these cells
    links_norm = [links_norm; n_user, nearest_macro; nearest_macro, n_user; n_user, nearest_pico; nearest_pico, n_user];
end

% Generate D2D users
% Start generating one of the devices
d2d_users_pos = [space(1)*rand(n_d2d, 1), space(2)*rand(n_d2d, 1)];
% Then generate another device near of the couple
for i=1:n_d2d
    d2d_users_pos(n_d2d+i,:) = [d2d_users_pos(i,1)+(dist_pico_x/5)*rand(1,1), d2d_users_pos(i,2)+(dist_pico_y/5)*rand(1,1)];
    n_device = i+n_pico+n_macro+n_sec;
    links_norm = [links_norm; n_device, n_device+n_d2d; n_device+n_d2d, n_device];
end
nodes_pos = [nodes_pos; d2d_users_pos];

% Calculate distance between every link in links_norm and links_fib
len_norm = size(links_norm,1);
len_fib = size(links_fib,1);
v_dist = zeros(len_norm+len_fib,1);
for i=1:len_fib
    v_dist(i) = sqrt((nodes_pos(links_fib(i,1),1)-nodes_pos(links_fib(i,2),1))^2+(nodes_pos(links_fib(i,1),2)-nodes_pos(links_fib(i,2),2))^2);
end
for i=1:len_norm
    v_dist(i+len_fib) = sqrt((nodes_pos(links_norm(i,1),1)-nodes_pos(links_norm(i,2),1))^2+(nodes_pos(links_norm(i,1),2)-nodes_pos(links_norm(i,2),2))^2);
end

% Create primary users
prim_pos = [space(1)*rand(n_prim, 1), space(2)*rand(n_prim, 1)];
links_prim = [];
v_dist_prim = [];

if n_prim ~= 0
    % Calculate distance between secondary and primary nodes
    dist_mat_prim = zeros(n_prim,2*n_d2d+n_sec);
    for i=1:n_prim
        for j=n_pico+n_macro+1:n_pico+n_macro+2*n_d2d+n_sec
            dist_mat_prim(i,j-n_pico-n_macro) = sqrt((nodes_pos(j,1)-prim_pos(i,1))^2+(nodes_pos(j,2)-prim_pos(i,2))^2);
        end
    end

    % Create primary links based on prim_links
    adj_mat_prim = zeros(n_prim,2*n_d2d+n_sec);
    for i=1:n_prim
        [dist_sorted, idx] = sort(dist_mat_prim(i,:));
        neighb = idx(1:prim_links);
        adj_mat_prim(i,neighb) = 1;
        for j=1:prim_links
            links_prim = [links_prim; i neighb(j)+n_pico+n_macro];
            v_dist_prim = [v_dist_prim; dist_sorted(j)];
        end
    end
end

% Create adjacency matrix for fiber links
n_bs = n_macro+n_pico;
adj_mat_fib = zeros(n_bs,n_bs);
for i=1:length(links_fib)
    adj_mat_fib(links_fib(i,1), links_fib(i,2)) = 1;
end
bs_pos = [macros_pos;picos_pos];

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
plot(macros_pos(:,1),macros_pos(:,2),'ro');
plot(picos_pos(:,1),picos_pos(:,2),'gp');
plot(sec_users_pos(:,1),sec_users_pos(:,2),'k*');
plot(d2d_users_pos(:,1),d2d_users_pos(:,2),'m^');
if n_prim ~= 0
    plot(prim_pos(:,1),prim_pos(:,2),'rx');
    tmp = 0;
    for i=1:n_prim
        for j=1:prim_links
            tmp = tmp+1;
            plot([prim_pos(i,1) nodes_pos(links_prim(tmp,2),1)], [prim_pos(i,2) nodes_pos(links_prim(tmp,2),2)], 'r:');
        end
    end
end
axis([0 space(1), 0 space(2)])
legend('Wireless links','Fiber links','Macrocells', 'Picocells', 'Secondary users', 'D2D users', 'Primary users', 'Interference links');
grid on;

end

