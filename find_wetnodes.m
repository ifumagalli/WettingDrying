function [wetnodes,littlewetnodes,frontnodes,frontwettednodes,firstdrynodes] = find_wetnodes(elements,u,c,g,tol,tol2,phase)

% [ wetnodes,frontnodes ] = FIND_WETNODES(elements,u,c,g,tol,phase)
% finds wet and front nodes and returns vectors with their indexes
% u = x-velocity
% c = celerity = 2*sqrt(g*wd)
% g = gravity acceleration
% phase =   'pred' 1st step of predictor-corrector algorithm
%           'corr' 2nd step of predictor-corrector algorithm
% NB:	- A node is said to be wet if it's a vertice of a triangle that has 
% 		  wd > tol in at least 2 vertices
%       - A wet node is said to be little wet if it has wd < tol2
%		- A node is said to be on the front if it is a vertice of a partially
%		  wet triangle: exactly 2 vertices have wd > tol
%       - frontwettednodes: front nodes on which wd < tol

switch phase
    case 'pred',
        tol = 2*sqrt(g*tol/10); % tol was on height, and we have celerity
        tol2 = 2*sqrt(g*tol2/10);
    case 'corr',
        tol = 2*sqrt(g*tol);
        tol2 = 2*sqrt(g*tol2);
    otherwise,
end

wetnodes = zeros(1,size(elements,2)*3);
littlewetnodes = [];%zeros(1,size(elements,2)*3);
frontnodes = [];
frontwettednodes = [];
for i=1:size(elements,2)
    %%% this cycle could be made more efficient
    v1 = elements(1,i); v2 = elements(2,i); v3 = elements(3,i);
    if ((c(v1)>tol) + (c(v2)>tol) + (c(v3)>tol))>1
%        ((u(v1)>c(v1)/2) + (u(v2)>c(v2)/2) + (u(v3)>c(v3)/2))>1
%     if ((abs(c(v1))>tol) + (abs(c(v2))>tol) + (abs(c(v3))>tol))>1
        wetnodes(3*i-2) = v1;
        wetnodes(3*i-1) = v2;
        wetnodes(3*i) = v3;
        if ((c(v1)>tol) + (c(v2)>tol) + (c(v3)>tol))<3
			frontnodes = [frontnodes; elements(:,i)];
        end
        if c(v1)<=tol
            frontwettednodes = [frontwettednodes; v1];
        elseif c(v1)<tol2
            littlewetnodes = [littlewetnodes; v1];
        end
        if c(v2)<=tol
            frontwettednodes = [frontwettednodes; v2];
        elseif c(v2)<tol2
            littlewetnodes = [littlewetnodes; v2];
        end
        if c(v3)<=tol
            frontwettednodes = [frontwettednodes; v3];
        elseif c(v3)<tol2
            littlewetnodes = [littlewetnodes; v3];
        end
    end
end

wetnodes = setdiff(unique(wetnodes),[0]);
    % setdiff is to remove 0: it's not an index
frontnodes = unique(frontnodes);
frontwettednodes = unique(frontwettednodes);
littlewetnodes = unique(littlewetnodes);
 
old_wetnodes = wetnodes;
for i=1:5
    new_wetnodes = [];
    for i=1:size(elements,2)
        v1 = elements(1,i); v2 = elements(2,i); v3 = elements(3,i);
    %     if ~isempty(union(union(find(old_wetnodes==v1),find(old_wetnodes==v2)),find(old_wetnodes==v3)))
        if (~isempty(find(old_wetnodes==v1))+~isempty(find(old_wetnodes==v2))+~isempty(find(old_wetnodes==v3)))>1
            new_wetnodes = [new_wetnodes v1 v2 v3];
        end
    end
    old_wetnodes = unique(new_wetnodes);
end
firstdrynodes = setdiff(old_wetnodes,wetnodes);

end