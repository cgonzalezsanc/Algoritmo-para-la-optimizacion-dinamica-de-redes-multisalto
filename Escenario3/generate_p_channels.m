function [h] = generate_p_channels(Graph, t, it, Sigma, nodes, p_users, d_Graph)

Points = t/it;
h = zeros (Points, p_users, nodes);

for i = 1:p_users
    for j = 1:nodes
        if (Graph(i, j) == 1) 
            h(:, i,j) = create_wless_chan(t, it, Sigma, d_Graph(i,j));
        end;
    end
end