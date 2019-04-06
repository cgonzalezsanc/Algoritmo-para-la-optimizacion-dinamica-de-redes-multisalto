% Función que crea los Maximal Independent Set

function [ws, S] = Independent_Sets(link, links_fb, Graph)

Num_link = length (link);

w0 = zeros (Num_link + length(links_fb), Num_link + length(links_fb));
ws_p = zeros (Num_link + length(links_fb), Num_link + length(links_fb));
ws = zeros (Num_link + length(links_fb), Num_link + length(links_fb));

% Un enlace es independiente (ortogonal) de otro e_ij si:
%
%   - El enlace no debe tener como origen o destino el nodo i
%   - El enlace no debe tener como origen o destino el nodo j
%   - El enlace no debe tener como destino un nodo vecino a i
%   - El enlace no debe tener como origen un nodo vecino a j
%
% En w0 se indica si un enlace es independiente de otro. Es el grafo de
% conflicto.
for i = 1:Num_link
    for j = 1:Num_link
         % Comprobaciones:
         %  - Que no sea origen i
         %  - Que no sea destino i
         %  - Que no sea origen j
         %  - Que no sea destino j
         %  - Que no sea destino un nodo vecino a i
         %  - Que no sea origen un nodo vecino a j
         if (link(j,1) == link(i,1)) || (link(j,2) == link(i,1)) || ((link(j,1) == link(i,2)) || (link(j,2) == link(i,2))) || (Graph(link(j,2), link(i,1)) == 1) || (Graph(link(i,2), link(j,1)) == 1)
            w0(i,j) = 0;
        else
            w0(i, j) = 1;
        end;
    end
end


for l=1:Num_link
    ws_p(l,1) = l;
end

k = Num_link;

c = 1;
cont = 1;

while cont == 1
   
    cont = 0;
    remove = 0;
    
    maximum = 0;
    while ws_p (c,maximum+1) > 0
        maximum = maximum + 1;
    end;
    
    for j = ws_p(c,maximum):Num_link
        
        p = 1;
        selected = 1;
        while ws_p (c,p) > 0
            if w0( ws_p(c,p), j) == 0
                selected = 0;
                break;
            end;
            p = p + 1;
        end;
        
        if selected == 1
            k = k+1;
            p = 1;
            while ws_p (c,p) > 0
                ws_p(k, p) = ws_p (c,p);
                p = p +1;
            end;
            ws_p(k,p) = j;
            cont = 1;
            remove = 1;
        end;
                
        if c < k
            cont = 1;
        end;
    end;
    
    if remove == 1
        ws_p(c,1) = -1;
    end;
    
    c = c + 1;
end;

s = 1;

for i = 1:k
    if ws_p(i,1) > -1
        ws(s,:) = ws_p(i,:);
        s = s + 1;
    end;
end;

S= 0;
while ws(S+1,1) > 0;
    S = S + 1;
    if S == length(ws)
        break;
    end;
end;

if ~isempty(links_fb)
    ws(ws~=0) = ws(ws~=0)+length(links_fb);
    for i=1:S
        a = 1;
        while ws(i,a) > 0
            a = a + 1;
        end;

        for z = 1:length(links_fb)
            ws(i,a + z - 1) = z;
        end;
    end;
end;

end
