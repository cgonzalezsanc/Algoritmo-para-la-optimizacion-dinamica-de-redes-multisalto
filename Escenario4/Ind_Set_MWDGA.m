function [mis, n_mis] = Ind_Set_MWDGA(m_adj, phi, links_norm, links_fib)
% Creates a Maximal Weighted Independent Set (MIS) using the Minimum 
% Weighted Degree Greedy Algorithm.
%
% Inputs:
%   - m_adj: Adjacency matrix
%   - phi: quality channel indicator
%
% Outputs:
%   - mis: MIS
%   - n_mis: Number of MIS

% Initializations
mis = [];
links = links_norm;
v_ind_link = 1:size(m_adj,1); % Vector which contains node indexes

% Exclude fiber links (these will be included in the MIS anyway)
for i=1:size(links_fib,1)
    orig = links_fib(i,1);
    dest = links_fib(i,2);
    phi(orig,dest) = 0;
end

% Start iterations
for i=1:size(m_adj,1)
    % Get best link to transmit
    bestQual = min(phi(:));
    if bestQual == 0
        break;
    end
    [orig,dest] = find(phi==bestQual);
    % Pick one if there are multiple best links
    if length(orig) > 1
       k = randi(length(orig));
       orig = orig(k);
       dest = dest(k);
    end
    % Look for that node in the conlfict graph
    link = find(links(:,1)==orig & links(:,2)==dest);
    % Add that link to the MIS
    mis(i) = v_ind_link(link);
    % Remove that link and its neighbours from the conflict graph
    v_neigh = find(m_adj(link,:) == 1);
    rm_links = unique([link, v_neigh]);
    m_adj(rm_links, :) = [];
    m_adj(:, rm_links) = [];
    v_ind_link(rm_links) = [];
    % Remove that link and its neighbours from the quality matrix (quality
    % equal to 0)
    for k = 1:length(rm_links)
        x1 = links(rm_links(k),1);
        x2 = links(rm_links(k),2);
        phi(x1,x2) = 0;
    end
    % Remove that link and its neighbours from the link vector
    links(rm_links, :) = [];
end

% Include fiber links
if ~isempty(links_fib)
    v_links_fib = 1:length(links_fib);
    mis = mis + length(links_fib);
    mis = [mis repmat(v_links_fib, size(mis,1), 1)];
    if isempty(mis)
        mis = v_links_fib;
    end
end

n_mis = size(mis,1);


end

