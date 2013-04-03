function frontnodes = find_front(elements,chi_wet)

% Finds the advancing shoreline nodes, given the characteristic function of
% wet nodes: a node is said to belong to the shoreline if it is a vertice
% of a triangle which has at least 1 wet vertices and not more than 2
% NB: DA CAMBIARE SE VOGLIAMO PRENDERE ANCHE LA FILA DI CELLE INTERNA
% 
% chi_wet:  chi_wet(i)=1 if i is a wet node
% front:    vector containing the index of shoreline nodes

frontnodes = [];
for i=1:max(size(elements))
    if (chi_wet(elements(1,i)) + chi_wet(elements(2,i)) + chi_wet(elements(3,i))>=1) ...
        & (chi_wet(elements(1,i))+ chi_wet(elements(2,i)) + chi_wet(elements(3,i))<=2)
        frontnodes = [frontnodes; elements(:,i)];
%         frontnodes(3*i-2) = elements(1,i);
%         frontnodes(3*i-1) = elements(2,i);
%         frontnodes(3*i) = elements(3,i); 
    end
    frontnodes = unique(frontnodes);
end

end