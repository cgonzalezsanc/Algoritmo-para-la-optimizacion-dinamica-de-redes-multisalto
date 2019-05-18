function [mis, n_mis] = Ind_Set_MaxProd(m_adj, phi, links_sets, links_fib)
% Creates a Maximal Weighted Independent Set (MWIS) using the min-sum
% algorithm (modification of Max Product algorithm)

% Declarations
m_ij = zeros(size(m_adj,1), size(m_adj,2), 2);      % messages at t
y_ij = zeros(size(m_adj,1), size(m_adj,2));         % gamma at t
y_ij_nw = zeros(size(m_adj,1), size(m_adj,2));      % gamma at t+1
x_i = zeros(size(m_adj,1),1);                       % scheduling variable  
max_t = 140;                                         % maximum number of iterations

% 0) Definitions (t=0)
m_ij = ones(size(m_ij));
y_ij = log(m_ij(:,:,1)./m_ij(:,:,2));          % must be zeros

% 1) Iterative algorithm
for t=1:max_t
    for i=1:size(m_adj,1)
        % Get link information
        link = links_sets(i,:);
        % Calculate params
        phi_ij = get_link_phi(link, phi);
        neighbors = get_neighbors(i, m_adj);
        for j=neighbors
            rest_neighbors = exclude_dest(neighbors,j);  % exclude j
            sum_y = sum_neighbors(rest_neighbors, y_ij(:,i));
            % Calculate gamma
            y_ij_nw(i,j) = max(phi_ij-sum_y,0);
        end
        % Compare phi and gamma
        if phi_ij > sum_neighbors(neighbors, y_ij(:,i))
            x_i(i) = 1;
        elseif phi_ij < sum_neighbors(neighbors, y_ij(:,i))
            x_i(i) = 0;
        else
            x_i(i) = 1/2;
        end
    end
    y_ij = y_ij_nw;
    y_ij_nw = zeros(size(y_ij_nw));
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

function phi_ij = get_link_phi(link, phi)
    phi_ij = -1*phi(link(1),link(2));
end

function sum_y = sum_neighbors(neighbors, y_ij)
    sum_y=0;
    for k=neighbors
        sum_y = sum_y+y_ij(k);
    end
end