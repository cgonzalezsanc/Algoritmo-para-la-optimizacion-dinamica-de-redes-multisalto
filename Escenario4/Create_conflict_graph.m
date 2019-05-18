function m_adj = Create_conflict_graph(links, Graph)

Num_link = length (links);
% Initialization of the adjacency matrix
m_adj = zeros(Num_link, Num_link);

for i = 1:Num_link
    for j = 1:Num_link
         % Check:
         %  - i can't be the origin
         %  - i can't be the destination
         %  - j can't be the origin
         %  - j can't be the destination
         %  - neighbors of i can't be the destination
         %  - neighbors of j can't be the origin
         if (links(j,1) == links(i,1)) || (links(j,2) == links(i,1)) || ((links(j,1) == links(i,2)) || (links(j,2) == links(i,2))) || (Graph(links(j,2), links(i,1)) == 1) || (Graph(links(i,2), links(j,1)) == 1)
            m_adj(i,j) = 1;
        else
            m_adj(i, j) = 0;
        end;
    end
end