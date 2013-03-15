% n~ con ~=nome di un vettore : lunghezza di quel vettore
%
% dof_v = gradi di libertà in velocità (sia u sia v) (indice vertice)
% dof_c = gradi di libertà in celerità (indice vertice)
% dof = [dof_v , ndof_v+dof_v, 2*ndof_v+dof_c]
% Dnodes = nodi Dirichlet (notazione come dof)
% wetnodes = nodi bagnati
% wetdof_uv = intersect(wetnodes,dof_uv_tot)
% wetdof_c = analogo
% wetdof = [wetdof_uv,wetdof_uv+ndof_v,wetdof_c+2*ndof_v]';
%
% NOTE:
% Si suppone che i nodi di Dirichlet siano gli stessi per u e v
%           (l'ha fatto già lui)