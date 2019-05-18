function [mis, n_mis] = Ind_Set_CMPA(m_adj, phi, links_sets, links_fib, eps, d)
% Creates a Maximal Weighted Independent Set (MWIS) using Convergent
% Message Passing Algorithm (CMPA).

% Initializations
n_nodes = size(m_adj,1);                % number of nodes
conf_graph_links = find(m_adj(:)==1);   % conflict graph links
n_links = length(conf_graph_links);     % number of links
y_ij = zeros(n_nodes, n_nodes);         % gamma at t
y_ij_nw = zeros(n_nodes, n_nodes);      % gamma at t+1
l_ij = zeros(n_nodes, n_nodes);         % lambda at t
l_ij_nw = zeros(n_nodes, n_nodes);      % lambda at t+1
w = zeros(n_nodes,1);                   % link weights
x_i = zeros(n_nodes,1);                 % scheduling variable

% Link weights
[r,c] = find(phi~=0);
for i = 1:length(r)
    % Get weight
    w_i = -1*phi(r(i),c(i));
    % Get link index
    link_idx = find(links_sets(:,1)==r(i) & links_sets(:,2)==c(i));
    % Assign this weight to weigth matrix
    w(link_idx) = w_i;
end

% 1) Calculation of lambda
% Initialization for t=0
for i=1:n_nodes
    for j=1:n_nodes
        l_ij(i,j) = max(w(i),w(j));
    end
end

max_t = 1*n_links;      % max number of iterations
for t = 1:3*max_t
    if mod(t,max_t)==0
        zzz=1;
    end
    % Get link index
    idx = rem(t-(max_t*floor((t-1)/max_t)),max_t+1);
    link_idx = conf_graph_links(idx);
    % Get origin (row) and destination (col)
    orig = floor(link_idx/n_nodes)+1;
    dest = rem(link_idx,n_nodes);
    if dest == 0
        orig = orig-1;
        dest = n_nodes;
    end
    % Calculate origin gamma
    neighbors = get_neighbors(orig, m_adj); % get origin neighbors
    rest_neighbors = exclude_dest(neighbors,dest);
    sum_y = sum_neighbors(rest_neighbors, y_ij(:,orig));
    y_ij_nw(orig,dest) = max(w(orig)-sum_y,0);
    % Calculate destination gamma
    neighbors = get_neighbors(dest, m_adj); % get destination neighbors
    rest_neighbors = exclude_dest(neighbors,orig);
    sum_y = sum_neighbors(rest_neighbors, y_ij(:,dest));
    y_ij_nw(dest,orig) = max(w(dest)-sum_y,0);
    % Calculate lambda
    a = y_ij_nw(orig,dest);
    b = y_ij_nw(dest,orig);
    l_ij(orig,dest) = (a+b+2*eps+sqrt((a-b)^2+4*(eps)^2))/2;
    if rem(t,max_t)==0
        y_ij = y_ij_nw;
    end
end

% 2) Calculation of MWIS
% colours are: 0 - gray, 1 - green, 2 - orange, 3 - red
colour = ones(n_nodes,1); % nodes are green initially

% Take a minimum node
[~,node] = min(w);
neighbors = get_neighbors(node, m_adj); % get node neighbors
sum_l = sum_neighbors(neighbors, l_ij(node,:));
if sum_l > w(node)+d
   colour(node) = 0;
   x_i(node) = 0;
end

% for i = 1:n_nodes
%     neighbors = get_neighbors(i, m_adj); % get node neighbors
%     sum_l = sum_neighbors(neighbors, l_ij(i,:));
%     if sum_l > w(i)
%        colour(i) = 0;
%        x(i) = 0;
%     end
% end

for ii = 1:3*n_nodes
    if mod(i,n_nodes) == 0
        i = 1;
    else
        i = i+1;
    end   
    if colour(i) == 1
        neighbors = get_neighbors(i, m_adj);
        neigh_gray = find(colour(neighbors)==0);
        for n = 1:length(neigh_gray)
            j=neigh_gray(n);
            if l_ij(i,j)>d
                x_i(i) = 1;
                colour(i) = 2;
                break;
            end
        end
        if any(colour(neighbors)==2)
            x_i(i)=0;
            colour(i)=0;
        end
    end
end

for i = find(colour==2)
   colour(i)=3;
   x_i(i) = 1;
end

mis = find(x_i==1).';
n_mis = 1;

end


function neighbors = get_neighbors(node, m_adj)
    neighbors = find(m_adj(node,:) == 1);
    neighbors(neighbors==node) = [];
end

function new_neighbors = exclude_dest(neighbors, dest)
    dest_ind = neighbors==dest;
    new_neighbors = neighbors;
    new_neighbors(dest_ind) = [];
end

function sum_a = sum_neighbors(neighbors, a_ij)
    sum_a=0;
    for k=neighbors
        sum_a = sum_a+a_ij(k);
    end
end