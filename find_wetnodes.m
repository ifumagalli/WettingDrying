function [ wetnodes ] = find_wetnodes(boundaries,c,g,tol,phase)
% [ wetnodes ] = FIND_WETNODES(boundaries&,c&,tol)
% finds the wet nodes and returns a vector with their indexes
% C = celerity 2*sqrt(g*wd)
% g = gravity acceleration
% phase =   'pred' 1st step of predictor-corrector algorithm
%           'corr' 2nd step of predictor-corrector algorithm

switch phase
    case 'pred',
        tol = 2*sqrt(g*tol); % tol was on height, and we have celerity
    case 'corr',
        tol = 2*sqrt(g*tol);
    otherwise,
end

% Wetting elements whose wet vertices are more than 2
% (only one step of front progression)
wetnodes = zeros(1,size(boundaries,2)*2);
for i=1:size(boundaries,2)
    %%% this cycle could be made more efficient
    v1 = boundaries(1,i); v2 = boundaries(2,i);
    if c(v1)>tol | c(v2)>tol
        wetnodes(2*i-1) = v1;
        wetnodes(2*i) = v2;
    end
end
    
wetnodes = unique(wetnodes);

end

