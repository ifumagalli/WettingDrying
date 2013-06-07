function [A,rhs] = assem_mat_vect_gio_stab_adapt(vertices,elements,boundaries,g,dt,un,vn,cn,una,vna,cna,theta,coeff,h0,h_k,f1,f2,f3,wdn,omega,time...
    ,un_old,vn_old,cn_old,nx_nodes,ny_nodes,frontnodes,frontwettednodes,littlewetnodes,drynodes,save_path)

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
[Gx2,Gy2] = ass_grad(vertices,elements,omega,1,2);

GxMom=Gx; GyMom=Gy;
% Per spegnere il termine con c/2 grad c nell'eq del momento sul fronte
% GxMom(frontnodes,:) = Gx2(frontnodes,:);
% GyMom(frontnodes,:) = Gy2(frontnodes,:);
% GxMom(littlewetnodes,:) = Gx2(littlewetnodes,:);
% GyMom(littlewetnodes,:) = Gy2(littlewetnodes,:);
% GxMom(drynodes,:) = Gx2(drynodes,:);
% GyMom(drynodes,:) = Gy2(drynodes,:);

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

% % Per spegnare la gravità sui frontnodes
% rhs_x(frontnodes) = 0*rhs_x(frontnodes);
% rhs_y(frontnodes) = 0*rhs_y(frontnodes);
% rhs_3(frontnodes) = 0*rhs_3(frontnodes);
% 
% % Per spegnare la gravità sui littlewetnodes
% rhs_x(littlewetnodes) = 0*rhs_x(littlewetnodes);
% rhs_y(littlewetnodes) = 0*rhs_y(littlewetnodes);
% rhs_3(littlewetnodes) = 0*rhs_3(littlewetnodes);

% Per la correzione di continuità di cui in [Horritt, eq (12)]
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

rhs = (M3 - (1-theta)*dt*A)*[un;vn;cn] + delta_s*(M_Lh - dt*(1-theta_s)*Lh)*[un;vn;cn] +dt*([rhs_x;rhs_y;rhs_3] - [rhs_fr;zeros(nov,1)] + delta_s*f_Lh) + CC_rhs;

% % Per modulare la stabilizzazione sul fronte
% delta_front = 0.5;
% Lh(frontnodes,:) = delta_front*Lh(frontnodes,:);
% Lh(:,frontnodes) = delta_front*Lh(:,frontnodes);
% Lh(frontnodes,frontnodes) = Lh(frontnodes,frontnodes) / delta_front;
% M_Lh(frontnodes,:) = delta_front*M_Lh(frontnodes,:);
% M_Lh(:,frontnodes) = delta_front*M_Lh(:,frontnodes);
% M_Lh(frontnodes,frontnodes) = M_Lh(frontnodes,frontnodes) / delta_front;
% f_Lh(frontnodes) = delta_front*f_Lh(frontnodes);

% matrice
A = M3 + theta*dt*A + delta_s*(M_Lh + dt*theta_s*Lh) + CC_A;



