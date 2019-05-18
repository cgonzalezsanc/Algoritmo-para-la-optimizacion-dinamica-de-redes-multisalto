function [link] = Get_links(Graph)

link = zeros(1,2);

pointer = 0;
for i = 1:length(Graph)
    for j = 1:length(Graph)
        if (Graph(i, j) == 1) 
            pointer = pointer + 1;
            link(pointer,:) = [i, j];
        end;
    end
end