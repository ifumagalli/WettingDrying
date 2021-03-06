function [A,rhs] = assem_mat_vect_gio_stab_adapt(vertices,elements,boundaries,g,dt,un,vn,cn,una,vna,cna,theta,coeff,h0,h_k,f1,f2,f3,wdn,omega,time...
    ,un_old,vn_old,cn_old,cn_vecchia,nx_nodes,ny_nodes,wetnodes,frontnodes,frontwettednodes,littlewetnodes,drynodes,save_path)

%
% assembla la matrice globale e il termine noto del sistema 
% con operazioni vettoriali

nov = size(vertices,2);
noe=size(elements,2);
% matrice di advezione: \int_T (\alpha U \cdot grad \phi_j) \phi_i
% A_adve_alpha = ass_adve(vertices,elements,un,vn,alpha,1);

% matrice di advezione: \int_T (U \cdot grad \phi_j) \phi_i
A_adve = ass_adve(vertices,elements,una,vna,omega,1,2);

% A_adve_delta = ass_adve(vertices,elements,un,vn,delta,2);

% matrici: 
% Gx = \int_T g d(\phi_j)/dx \phi_i
% Gy = \int_T g d(\phi_j)/dy \phi_i 
[Gx,Gy] = ass_grad(vertices,elements,omega,cna*0.5,2);

GxMom=Gx; GyMom=Gy;
% Per spegnere il termine con c/2 grad c nell'eq del momento sul fronte
% GxMom(frontnodes,:) = 0;
% GyMom(frontnodes,:) = 0;
% GxMom(littlewetnodes,:) = 0;
% GyMom(littlewetnodes,:) = 0;
% GxMom(drynodes,:) = 0;
% GyMom(drynodes,:) = 0;

% matrici: 
% Hx = \int_T H d(\phi_j)/dx \phi_i
% Hy = \int_T H d(\phi_j)/dy \phi_i 
% [Hx,Hy] = ass_grad(vertices,elements,1,cn/2,1);

% matrice di massa

M=ass_massa(vertices,elements,omega,1,2);

% matrice associata alle resistenze al fondo

% Fr = ass_massa(vertices,elements,omega,coeff,2);

% 
% matrice nulla
Zero = sparse(nov,nov);

% rhs_fr = [Fr  , Zero  ;
%          Zero , Fr     ]* [una;vna];
% rhs_fr = [coeff'.*una;coeff'.*vna];

[rhs_x,rhs_y,rhs_3] = ass_rhs_g(vertices,elements,f1,f2,f3,4);
[rhs_frx,rhs_fry] = ass_rhs_fr(vertices,elements,coeff,un,vn,4);
rhs_fr = [rhs_frx;rhs_fry];
% [Fx,Fy] = ass_grad(vertices,elements,1,c0*0.5,1);
% rhs_x = Fx*c0;
% rhs_y = Fy*c0;

% contorni impermeabili

[N11,N12,N21,N22] = ass_normali_penalty(vertices,elements,nx_nodes,ny_nodes,2);
gamma = 1e10;

N = gamma*[N11 , N12 , Zero;
           N21 , N22 , Zero;
           Zero, Zero, Zero];

 
% matrice di massa a blocchi
M3 = [M   , Zero, Zero ;
     Zero , M   , Zero ;
     Zero , Zero, M];


% matrice finale
% Dx=Hx-F_x;
% Dy=Hy-F_y;

A =  [A_adve     ,Zero     ,GxMom  ;
      Zero       ,A_adve   ,GyMom  ;
      Gx         ,Gy       ,A_adve];

% A = A + [Fr         , Zero      ,   Zero  ;
%          Zero       , Fr        ,  Zero  ;
%          Zero       , Zero      ,  Zero  ];
%      
% 
A = A + N;

%parametro e matrici stabilizazzione
% Porre delta_s=0 per spegnere tutta la stabilizzazione: sia streamline
% diffusion, sia shock capturing (per spegnere solo shock capturing agire
% in assem_mi, per spegnere solo streamline diffusion BOH)
% delta_s = max(h_k);
delta_s = 1;
% delta_s = 0;
theta_s = theta;
[M_Lh,Lh,f_Lh] = assem_mat_GLS_simm1(vertices,elements,boundaries,g,una,vna,cna,coeff,f1,f2,f3,theta,dt,h_k,omega,time,un,vn,cn,un,vn,drynodes,save_path);

% % Per spegnare la gravit� sui frontnodes
% rhs_x(frontnodes) = 0*rhs_x(frontnodes);
% rhs_y(frontnodes) = 0*rhs_y(frontnodes);
% rhs_3(frontnodes) = 0*rhs_3(frontnodes);
% 
% % Per spegnare la gravit� sui littlewetnodes
% rhs_x(littlewetnodes) = 0*rhs_x(littlewetnodes);
% rhs_y(littlewetnodes) = 0*rhs_y(littlewetnodes);
% rhs_3(littlewetnodes) = 0*rhs_3(littlewetnodes);

% Per la correzione di continuit� di cui in [Horritt, eq (12)]
%   sui nodi completamente bagnati degli elementi parzialmente bagnati
%   aggiungere eta*A_elem/2*(dc1/dt+dc2/dt)
eta = 0;
area = nanmean(pdetrg(vertices,elements));
temp = setdiff(frontnodes,frontwettednodes);
[temp1,temp2] = meshgrid(frontnodes,frontnodes);
CC_A = eta*area/2*M.*sparse(temp1,temp2,1,nov,nov);
CC_A = diag(CC_A*ones(nov,1)); % lumping
CC_A = [Zero Zero Zero; Zero Zero Zero; Zero Zero CC_A];
CC_rhs = CC_A*[0.*un;0.*vn;cn];

% % Per modulare la stabilizzazione sul fronte
% delta_front = 0.5;
% Lh(frontnodes,:) = delta_front*Lh(frontnodes,:);
% Lh(:,frontnodes) = delta_front*Lh(:,frontnodes);
% Lh(frontnodes,frontnodes) = Lh(frontnodes,frontnodes) / delta_front;
% M_Lh(frontnodes,:) = delta_front*M_Lh(frontnodes,:);
% M_Lh(:,frontnodes) = delta_front*M_Lh(:,frontnodes);
% M_Lh(frontnodes,frontnodes) = M_Lh(frontnodes,frontnodes) / delta_front;
% f_Lh(frontnodes) = delta_front*f_Lh(frontnodes);

% Per imporre u=umax sui dry
% ESITO PROVE
% 0-1-2-3-4-5 : u oscilla con valori anche negativi. Per forza, non ho
% messo il - davanti al calcolo di val!
% 0bis: u oscilla di bestia
% 1bis-3bis-4bis: u esplode
% 2bis: u oscilla con valori anche negativi
% 5bis-5bis_ogni2Corr: u esplode, per� fino a 0.8 h va bene, e poi il problema si genera da y=0
% 5bis_ogni3Corr: come 5 bis, ma con problema che parte anche da y=max, ma un pochino pi� tardi
% 5bis_ogni4Corr: esplode->salvata lo stesso
% 5bis_ogni5Corr: OK->salvata
[dx_cn_vecchia,dy_cn_vecchia]=pdegrad(vertices,elements,cn_vecchia);
dx_cn_vecchia = pdeprtni(vertices,elements,dx_cn_vecchia);
dy_cn_vecchia = pdeprtni(vertices,elements,dy_cn_vecchia);
val = zeros(nov,1);
% temp = setdiff(wetnodes,frontnodes); % prova 2-2bis
temp = setdiff(frontnodes,frontwettednodes); % prova 3-4 3bis
% val(frontnodes) = (cna(frontnodes)-cn_vecchia(frontnodes))./(dt*dx_cn_vecchia(frontnodes)+eps); %prove <=1
% val(frontnodes) = -(cna(frontnodes)-cn_vecchia(frontnodes))./(dt*dx_cn_vecchia(frontnodes)+eps); %prove 0bis-1bis
% val(temp) = (cna(temp)-cn_vecchia(temp))./(dt*dx_cn_vecchia(temp)+eps); %prova 2
val(temp) = -(cna(temp)-cn_vecchia(temp))./(dt*dx_cn_vecchia(temp)+eps); %prova 2bis
ys = unique(vertices(2,:));
for i=1:length(ys)
%     [~,idx] = min(abs(vertices(2,:)-(ymax+ymin)/2));
%     ymed = vertices(2,idx);
%     idxs_ymed = find(vertices(2,:) == ymed);
    idxs = find(vertices(2,:) == ys(i));
%     idx_front = max(intersect(idxs,frontnodes)); % prova 0-0bis
%     idx_front = min(intersect(idxs,frontnodes)); % prova 1-1bis
    idx_front = max(intersect(idxs,temp)); % prova 2-3-5-2bis-5bis
%     idx_front = min(intersect(idxs,temp)); % prova 4-4bis
    val(idxs(idxs >= idx_front)) = val(idx_front);
end
val'
% val = max(un(setdiff(wetnodes,frontnodes)));
% [V11,V12,V21,V22] = u_fixval_penalty(vertices,elements,val,2);
V11 = M;
V12 = 0.*V11;
V21 = 0.*V11;
V22 = 0.*V11;
% V(un(setdiff([1:nov],drynodes)),un(setdiff([1:nov],drynodes)))=0*V(un(setdiff([1:nov],drynodes)),un(setdiff([1:nov],drynodes)));
% V(un(setdiff([1:nov],drynodes)),:)=0*V(un(setdiff([1:nov],drynodes)),:);
% V(:,un(setdiff([1:nov],drynodes)))=0*V(:,un(setdiff([1:nov],drynodes)));
% V11(setdiff(wetnodes,frontnodes),:)=0; % prova <=4, (<=4)bis
V11(wetnodes,:)=0; % prova 5-5bis
if mod(floor(time/dt),4)~=0 
    V11 = Zero;             
end                         
% V(:,wetnodes)=0*V(:,wetnodes);
V = gamma*[V11 , V12 , Zero;
           V21 , V22 , Zero;
           Zero, Zero, Zero];
% rhs_V = V*[val*ones(nov,1); zeros(nov,1); zeros(nov,1)];
rhs_V = V*[val; zeros(nov,1); zeros(nov,1)];
rhs_V(setdiff(wetnodes,frontnodes)) = 0;

% termine noto e matrice IN QUESTO ORDINE
rhs = (M3 - (1-theta)*dt*A)*[un;vn;cn] + delta_s*(M_Lh - dt*(1-theta_s)*Lh)*[un;vn;cn] +dt*([rhs_x;rhs_y;rhs_3] - [rhs_fr;zeros(nov,1)] + delta_s*f_Lh) + CC_rhs + rhs_V;
A = M3 + theta*dt*A + delta_s*(M_Lh + dt*theta_s*Lh) + CC_A + V;

