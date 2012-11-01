vediamo ora
ciao by ele
ciao moreis sono l'ivan

% Contenuto della mesh:
% vertices ->   x
%               y
%               indice sottodominio
% boundaries -> indice vertice 1 (suo numero di colonna in vertices)
%               indice vertice 2
%               altro x4
% elements ->   indice vertice 1
%               indice vertice 2
%               indice vertice 3
% NB il tutto viene riordinato secondo convenienza (percorrenza antioraria)
% Noi useremo solo vertices(1:2,:),boundaries(1:2,:),elements(:,:)