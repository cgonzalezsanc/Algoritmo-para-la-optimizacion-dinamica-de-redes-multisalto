function [h] = generate_channels(Graph, t, it, Sigma, d_Graph)

Points = t/it;
h = zeros(Points, length(Graph), length(Graph));

for i = 1:length(Graph)
    for j = 1:length(Graph)
        if (Graph(i, j) == 1) 
            % Wireless link
            h(:, i,j) = create_wless_chan(t, it, Sigma, d_Graph(i,j));
        elseif Graph(i, j) == 2
            % Fiber link
            h(:, i,j) = Sigma;
        end;
    end
end