function [ wetnodes ] = find_wetnodes(elements,c,g,tol,phase)
% [ wetnodes ] = FIND_WETNODES(elements,c,tol)
% finds the wet nodes and returns a vector with their indexes
% C = celerity 2*sqrt(g*wd)
% g = gravity acceleration
% phase =   'pred' 1st step of predictor-corrector algorithm
%           'corr' 2nd step of predictor-corrector algorithm

switch phase
    case 'pred',
        tol = 2*sqrt(g*tol/10); % tol was on height, and we have celerity
    case 'corr',
        tol = 2*sqrt(g*tol);
    otherwise,
end

% Wetting elements whose wet vertices are at least 2
% (only one step of front progression)
wetnodes = zeros(1,size(elements,2)*3);
for i=1:size(elements,2)
    %%% this cycle could be made more efficient
    v1 = elements(1,i); v2 = elements(2,i); v3 = elements(3,i);
    if ((c(v1)>tol) + (c(v2)>tol) + (c(v3)>tol))>1
        wetnodes(3*i-2) = v1;
        wetnodes(3*i-1) = v2;
        wetnodes(3*i) = v3;
    end
end
    
wetnodes = setdiff(unique(wetnodes),[0]);
    % setdiff is to remove 0: it's not an index

end

