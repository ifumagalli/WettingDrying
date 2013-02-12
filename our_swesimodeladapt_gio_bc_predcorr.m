 

function our_swesimodeladapt_gio_bc_predcorr(vertices,elements,boundaries,tspan,tstart,theta)

% our_swesimodeladapt_gio_bc_predcorr(vertices,elements,boundaries,tspan,ts
% tart,theta)


%@>% = inizio modifica Bulgarello-Fumagalli
%@<% = fine modifica Bulgarello-Fumagalli

% h0 = bottom elevation
% wd = water depth
% h = h0+wd

%ELENA%
%save_path = 'C:\Users\Elenucci@\Dropbox\Elena&Ivan\AnEDP2\Progetto ANEDP2\saves\'
%IVAN%
save_path = 'C:\Users\Ivan\Dropbox\Elena&Ivan\AnEDP2\Progetto ANEDP2\saves\'

coarse=3; % alpha sempre =1 %??? unused?

dt = (tspan(2)-tspan(1))/tspan(3)
[h_k]=dim(vertices,elements); % diam(cerchio circoscritto)
% pdesurf(vertices,elements,h_k)
ho_k = h_k;
h_k_m=max(h_k)
% threshold nodi wet
wdtol = 1e-1;
deltat_deltax = dt/h_k_m 
[ar]=pdetrg(vertices,elements); % aree

% Dimensions
[d,noe]   = size(elements);
[d,nov]   = size(vertices);
[d,nside] = size(boundaries);

nbn = 2;

identita = eye(3*nov,3*nov);

[sk,lk,r1,r2]=stretching(vertices,elements);
% l2=lk;
% h_k = 2*lk;
pdesurf(vertices,elements,h_k)
% pause

% Common vectors
x = (vertices(1,:))';
y = (vertices(2,:))';
it1 = elements(1,:);
it2 = elements(2,:);
it3 = elements(3,:);


[xq,yq,wq,phiq] = basis_on_quad2D('P1',4);
xqq = phiq(1,:)'*x(it1)' + phiq(2,:)'*x(it2)' + phiq(3,:)'*x(it3)';
yqq = phiq(1,:)'*y(it1)' + phiq(2,:)'*y(it2)' + phiq(3,:)'*y(it3)';
% gravity acceleration, roughness coefficient, Coriolis coefficient
[g,n,l] = swecoeffgRl_gio; 

% Bathimetry: bottom elevation
% h0 = batswe2D_gio(x,y);
h0 = 1.61e-7*(x.^2 + y.^2);% parabolic bowl
b = abs(max(y)-min(y));

% Termine noto pendenza
[f_x,f_y]=pdegrad(vertices,elements,h0); %grad(h0) a centro cella


%f1 = -g*f_x; 
%f2 = -g*f_y;
%f3 = 0.*f1;
f1 = -g*f_x;
f1 = f1(ones(length(wq),1),:);
f2 = -g*f_y;
f2 = f2(ones(length(wq),1),:);
f3 = 0.*f1;

% % termine noto esatto caso bump
% h0qq = 0.2*exp(-0.05*(xqq-10).^2);
% dh0x = 0.02.*exp(-0.05.*(xqq-10).^2).*(10-xqq);
% 
% 
% dh0x = 2*wq*dh0x;
% 
% f1 = -g*dh0x;
% f2 = f1.*0;
% f3 = f1.*0;
% 
% pdesurf(vertices,elements,h0), pause
% pdesurf(vertices,elements,dh0x), pause

% termine noto esatto caso parabolic 


% rqq = sqrt(xqq.^2 + yqq.^2); 
% dh0x = 2*1.6*1e-7*xqq;
% f1 = -g*dh0x;
% dh0y = 2*1.6*1e-7*yqq;
% f2 = -g*dh0y;
% f3 = 0.*f1;

% termine noto esatto caso inventato hyppie
% 
% f1 = (xqq.^2 - xqq).*(2.*xqq-1) + 0.5*(10 + 2*cos(xqq)+tanh(10*(xqq - yqq)))...
%     .*(10*(1-tanh(10.*(xqq-yqq)).^2) - 2.*sin(xqq));
% 
% f2 = (yqq.^2 - yqq).*(2.*yqq-1) + 0.5*(10 + 2*cos(xqq)+tanh(10.*(xqq - yqq)))...
%      .*(-10.*(1-tanh(10.*(xqq-yqq)).^2));
% 
% f3 = (xqq.^2 - xqq).*(10.*(1-tanh(10.*(xqq-yqq)).^2) - 2*sin(xqq)) ...
%     +(yqq.^2 - yqq).*(-10.*(1-tanh(10.*(xqq-yqq)).^2))+ 0.5.*(10 + 2*cos(xqq)+tanh(10.*(xqq - yqq)))...
%     .*((2.*xqq-1) + (2.*yqq-1));


% f1 = (xqq.^2 - xqq).*(2.*xqq-1) + 0.5*(10 + 2*cos(xqq)+tanh(10*(xqq - yqq))).*(10*(1-tanh(10.*(xqq-yqq)).^2) - 2.*sin(xqq));
% 
% f2 =  0.5*(10 + 2*cos(xqq)+tanh(10.*(xqq - yqq))).*(-10.*(1-tanh(10.*(xqq-yqq)).^2));
% 
% f3 = (xqq.^2 - xqq).*(10.*(1-tanh(10.*(xqq-yqq)).^2) - 2*sin(xqq)) + 0.5.*(10 + 2*cos(xqq)+tanh(10.*(xqq - yqq)))...
%     .*(2.*xqq-1);
% f1 =  0.5*(10+2*cos(xqq)+tanh(10*(xqq - yqq))).*(10*(1-(tanh(10.*(xqq-yqq))).^2) - 2*sin(xqq));
% 
% f2 =  0.5*(10+2*cos(xqq)+tanh(10*(xqq - yqq))).*(-10.*(1-(tanh(10*(xqq-yqq))).^2));
% % 
% f3 = 0.*f2;
% f1 =  0.5*(10 + xqq.^2 + yqq.^2).*(2*xqq);
% 
% f2 =  0.5*(10 + xqq.^2 + yqq.^2).*(2*yqq);
% 
% f3 = 0.*f1;

% f1 = 2*wq*f1;
% f2 = 2*wq*f2;
% f3 = 2*wq*f3;

%griglia
% [sk,lk,r1,r2]=stretching(vertices,elements);
% l2=lk;
% l1=sk.*lk;

% Initial conditions

[un,vn,cin,wdin]=our_iniswe2D_gio(x,y,vertices,elements,tspan(1),g,h0);
wdtol = wdtol*max(wdin) % so that we have a relative tolerance

size(wdin)
size(h0)

hn = +wdin + h0; % altezza su piano di riferimento
wdn = wdin;
cn = cin;

wdn_el = pdeintrp(vertices,elements,wdn); % nodi -> centro cella
Vol(1)=sum(wdn_el.*ar)

pdesurf(vertices,elements,hn)
hold on
pdesurf(vertices,elements,h0)
hold off
figure(1001);
pdeplot(vertices,boundaries,elements,'xydata',wdn,'contour','on'), axis equal
title('Initial condition - water depth only')%strcat('t = ',num2str(t), ' - Corrector'))
pause;
% pdeplot(vertices,boundaries,elements,'zdata',hn)
% hold on
% pdesurf(vertices,elements,h0)
% hold off
% pause
   

% Normals
nx = [];
ny = [];
nx_edge = [];
ny_edge = [];

disp('ordina edge')
tic
% simike
[boundaries] = ordina_edge(boundaries,elements);

toc

disp('calcolo normali')
tic
for iside = 1:nside
    i1 = boundaries(1,iside);
    i2 = boundaries(2,iside);
    l = sqrt((vertices(1,i1)-vertices(1,i2))^2+(vertices(2,i1)-vertices(2,i2))^2);
    n1 = -(vertices(2,i1)-vertices(2,i2))/l;
    n2 = -(vertices(1,i2)-vertices(1,i1))/l;
    nx = [nx, n1 n1];
    ny = [ny, n2 n2];
    nx_edge = [nx_edge, n1];
    ny_edge = [ny_edge, n2];
end
toc

% vettore condizioni al contorno
% 1 è u.n=0, 2 è dirichlet o niente, problema ponte 

%  bc=[1,2,1,2];
% bc=[1,1,1,1,1,2,1,1,1,1,1,2];
% bc = [1,2,1,2]; 
% bc=[1,1,1,1,1,1,2]; %sudden
% bc=[1,1,2,1,2]; curva rett
% bc = [1,2,1,2,1,2];% conf
%  bc = [1,1,2,1,2];% obl
bc=[1,1,1,1]; %Radial Dam
% lati inflow ed outflow
in = 0;
out = 0;
in1 = 0;


disp('trova nodi Dir')
% tic
% Dirichlet boundaries for xi

DirDof_for_c0 = [find(boundaries(3,:) == out)];
DirDof_for_c1 = [find(boundaries(3,:) == 0)];%,find(boundaries(5,:)==4)];%Questo comando dovrà agire sulla riga 5, quella dove è identificato il lato
DirDof_for_c = [unique(boundaries(1:nbn,DirDof_for_c0));unique(boundaries(1:nbn,DirDof_for_c1))];% nodi di Dirichlet (prima eran dei lati, ora dei nodi)
%%% OSS: unique restituisce un vettore colonna

% Dirichlet boundaries for u and v
DirDof_for_uv_in0 = find(boundaries(3,:) == in);
DirDof_for_uv_in1 = find(boundaries(3,:) == in1);%Questo comando dovrà agire sulla riga 5, quella dove è identificato il lato
DirDof_for_uv_in = [unique(boundaries(1:nbn,DirDof_for_uv_in0));unique(boundaries(1:nbn,DirDof_for_uv_in1))];

% Boundary conditions
% contorni impenetrabili

% definisco le normali associate ai nodi...come quelle dei lati di bordo
% di cui i nodi sono primi nodi 
nx_nodes = sparse(nov,1);
ny_nodes = sparse(nov,1);
nx_nodes1 = sparse(boundaries(1,:),1,nx_edge,nov,1);
ny_nodes1 = sparse(boundaries(1,:),1,ny_edge,nov,1);
nx_nodes2  = sparse(boundaries(2,:),1,nx_edge,nov,1);
ny_nodes2  = sparse(boundaries(2,:),1,ny_edge,nov,1);

% definisco i vettori 
nx_nodes = nx_nodes + (nx_nodes1 + nx_nodes1)/2;
ny_nodes = ny_nodes + (ny_nodes1 + ny_nodes1)/2;
nx_nodes(DirDof_for_uv_in) = 0;
nx_nodes(DirDof_for_c) = 0;

ny_nodes(DirDof_for_uv_in) = 0;
ny_nodes(DirDof_for_c) = 0;

nodesxy =[];

% Sistemo alla vecchia maniera (codice precedente) il caso in cui u.n = 0
% sia imposto su tutti i lati
if sum(bc) == length(bc)
    for i=1:length(bc)  
        sides_u_dot_n = find(boundaries(7,:) == 1); % lati u_dot_n
    
isidesx = find( abs(nx_edge(sides_u_dot_n)) <   abs(ny_edge(sides_u_dot_n)) );
isidesy = find( abs(nx_edge(sides_u_dot_n)) >=  abs(ny_edge(sides_u_dot_n)) );

nodesx = union( boundaries(1,sides_u_dot_n(isidesx)),...
                boundaries(2,sides_u_dot_n(isidesx)) );
nodesy = union( boundaries(1,sides_u_dot_n(isidesy)),...
                boundaries(2,sides_u_dot_n(isidesy)) );
% nodi in comune
 nodesxy = intersect(nodesx,nodesy);
                                                                                
nnzy  = length(nodesy);
nnzx  = length(nodesx);
nnzxy = length(nodesxy);  

    end
else
%  Cerco di sistemare i nodi d'angolo che sono condivisi da lati u.n = 0 e Dir.
%  L'orientamento dei vettori normali va sempre verificato
    
for i=1:length(bc)  
        j=i;
        jp2=mod(i+length(bc)-1,length(bc)); 
        jp2=mod(i+2,length(bc))+1;% lato i-1
 % Se bc(i) = 1 e bc(i-1) = 0: quindi da wall a Dirichlet, metto la normale del nodo d'angolo 
 % come quella associata al secondo nodo       
 if bc(j)-bc(jp2) > 0
    nodes_j = [find(boundaries(3,:) == i)];
    nodes_j = [unique(boundaries(1:nbn,nodes_j))];
    nodes_jp2 = [find(boundaries(3,:) == jp2)];%,find(boundaries(5,:)==4)];%Questo comando dovrà agire sulla riga 5, quella dove è identificato il lato
    nodes_jp2 = [unique(boundaries(1:nbn,nodes_jp2))];
    corner = intersect(nodes_j,nodes_jp2);
    nx_nodes(corner) = nx_nodes2(corner);
    ny_nodes(corner) = ny_nodes2(corner);
 % DirDof_for_c = setdiff(DirDof_for_c ,corner);
 % Se bc(i) = 0 e bc(i-1) = 1: quindi da Dirichlet a wall, metto la normale del nodo d'angolo 
 % come quella associata al primo nodo
 elseif  bc(j)-bc(jp2) < 0   
    nodes_j = [find(boundaries(3,:) == i)];
    nodes_j = [unique(boundaries(1:nbn,nodes_j))];
    nodes_jp2 = [find(boundaries(3,:) == jp2)];%,find(boundaries(5,:)==4)];%Questo comando dovrà agire sulla riga 5, quella dove è identificato il lato
    nodes_jp2 = [unique(boundaries(1:nbn,nodes_jp2))];
    corner = intersect(nodes_j,nodes_jp2);
    nx_nodes(corner) = nx_nodes1(corner);
    ny_nodes(corner) = ny_nodes1(corner);
 %  DirDof_for_c = setdiff(DirDof_for_c ,corner);
 end
end
end


pdemesh(vertices,boundaries,elements), axis equal
hold on
quiver(vertices(1,:)',vertices(2,:)',nx_nodes,ny_nodes)
plot(vertices(1,DirDof_for_c),vertices(2,DirDof_for_c),'o')
plot(vertices(1,DirDof_for_uv_in),vertices(2,DirDof_for_uv_in),'x')
hold off
% % pause

dof_uv_in = setdiff([1:nov],DirDof_for_uv_in);%differenza insiemistica dei due vettori
ndof_uv_in = length(dof_uv_in);

dof_c = setdiff([1:nov],DirDof_for_c); %differenza insiemistica dei due vettori, nodi "interni"

% PER RISALTO
% dof_c = setdiff(dof_c,DirDof_for_uv_in);
dof_v = dof_uv_in;
ndof_v = length(dof_v);
ndof_c = length(dof_c);

% Dof and Dirichlet nodes
dof = [dof_v,dof_v+nov,dof_c+2*nov]; %veri dof
Dnodes = [DirDof_for_uv_in',DirDof_for_uv_in'+nov, DirDof_for_c'+2*nov];

% PER RISALTO
% Dnodes = [DirDof_for_uv_in',DirDof_for_uv_in'+nov, DirDof_for_c'+2*nov, DirDof_for_uv_in'+2*nov]; % dof di bordo

edges_in = find(boundaries(3,:) == in); 
edges_in1 = find(boundaries(3,:) == in1); 
edges_out = find(boundaries(3,:) == out);
nodes0 = unique(boundaries(1:nbn,edges_in));
nodes1 = unique(boundaries(1:nbn,edges_in1));
nodes2 = unique(boundaries(1:nbn,edges_out));

toc


% Inizializzo il calcolo delle portate
if sum(bc)~=length(bc)  % cioè se non è tutto muro (u.n=0) intorno
edges_in = find(boundaries(3,:) == in); 
nodes0 = unique(boundaries(1:nbn,edges_in));
dist_in = abs(vertices(2,boundaries(2,edges_in))-vertices(2,boundaries(1,edges_in)));
qx_in=0.5*(wdn(nodes0(1:end-1))+wdn(nodes0(2:end))).*dist_in'.*0.5.*(un(nodes0(1:end-1))+un(nodes0(2:end)));
Q_in(1)= sum(qx_in);

edges_in1 = find(boundaries(3,:) == in1); 
nodes1 = unique(boundaries(1:nbn,edges_in1));
dist_in1 = abs(vertices(1,boundaries(2,edges_in1))-vertices(1,boundaries(1,edges_in1)));
qx_in1=0.5*(wdn(nodes1(1:end-1))+wdn(nodes1(2:end))).*dist_in1'.*0.5.*(un(nodes1(1:end-1))+un(nodes1(2:end)));
qy_in1=0.5*(wdn(nodes1(1:end-1))+wdn(nodes1(2:end))).*dist_in1'.*0.5.*(vn(nodes1(1:end-1))+vn(nodes1(2:end)));

Q_in1(1)= sum(qy_in1);


edges_out = find(boundaries(3,:) == out); 
nodesout = unique(boundaries(1:nbn,edges_out));
dist_out = abs(vertices(2,boundaries(2,edges_out))-vertices(2,boundaries(1,edges_out)));
qx_out=0.5*(wdn(nodesout(1:end-1))+wdn(nodesout(2:end))).*dist_out'.*0.5.*(un(nodesout(1:end-1))+un(nodesout(2:end)));
qy_out=0.5*(wdn(nodesout(1:end-1))+wdn(nodesout(2:end))).*dist_out'.*0.5.*(vn(nodesout(1:end-1))+vn(nodesout(2:end)));
Q_out(1)= sum(qx_out);

end

% Inizializzo controllo di convergenza

[xq,yq,wq,phiq] = basis_on_quad2D('P1',4);
it1 = elements(1,:);
it2 = elements(2,:);
it3 = elements(3,:);
u_qq = phiq(1,:)'*un(it1)' + phiq(2,:)'*un(it2)' + phiq(3,:)'*un(it3)';
v_qq = phiq(1,:)'*vn(it1)' + phiq(2,:)'*vn(it2)' + phiq(3,:)'*vn(it3)';
h_qq = phiq(1,:)'*hn(it1)' + phiq(2,:)'*hn(it2)' + phiq(3,:)'*hn(it3)';

% Inizializzo variabili old e oold, che mi servono per stab e controllo di convergenza

unold = un;
vnold = vn;
cnold = cn;
unoold = un;
vnoold = vn;
cnoold = cn;

disp('Starting time loop')
% Time cicle
tspan(1),dt,tspan(2)
for t = tspan(1)+dt:dt:tspan(2)
    t
p=int64(t/dt);    
p1=int64(t/dt+1);

disp('assemblaggio')
tic 
    % matrice di stiffness e rhs 
  
    nod_zero=find(wdn==0);
    nod_frict= setdiff([1:nov],nod_zero);
    
    sigma_res(nod_frict) = n^2*g*sqrt(un(nod_frict,1).^2+vn(nod_frict,1).^2)./wdn(nod_frict,1).^(4/3);% coefficiente per resistenze al fondo
    sigma_res(nod_zero)= 0;
        % nella tesi sigma_res=S_f * g/u_vettore
        
%   sigma_res_t = pdeintrp(vertices,elements,sigma_res');  % ricavo un valore per ogni elemento

[aglo,rhs] = assem_mat_vect_gio_stab_i(vertices,elements,boundaries,g,dt,un,vn,cn,theta,sigma_res,h0,h_k,f1,f2,f3,wdn,1,t...
    ,unoold,vnoold,cnoold,nx_nodes,ny_nodes); 
figure(101), spy(aglo); title('matrice aglo')
if ~isempty(nodesxy)
%     aglo(nodesxy,:) = sparse(1:nnzxy,nodesxy,ones(1,nnzxy),nnzxy,3*nov,...
%                                  nnzxy);
%     aglo(nodesxy+nov,:) = sparse(1:nnzxy,nodesxy+nov,ones(1,nnzxy),...
%                                      nnzxy,3*nov,nnzxy);
%     rhs([nodesxy,nodesxy+nov]) = 0;   
end
toc
    disp('bc dir')
    tic
    % Dirichlet boundary conditions for xi
    [cDir0,cDir1,vDir] = xidirswe2D(vertices(1,DirDof_for_c),vertices(2,DirDof_for_c),boundaries,vertices,h0(DirDof_for_c),t,g,out,0);
    
    [uDir_in0,vDir_in0,uDir_in1,vDir_in1,cDir_in] = uvdir_inswe2D(vertices(1,DirDof_for_uv_in),vertices(2,DirDof_for_uv_in),t,wdn,DirDof_for_uv_in,g,in,in1,boundaries,vertices);
        
    disp('soluz lin sys')
%     tic
 
  
% Predictor
uDir_in = [uDir_in0,uDir_in1];

vDir_in = [vDir_in0,vDir_in1];
cDir = [cDir0,cDir1];

%@>%
disp('trova nodi wet')
tic
% Find wetnodes (da tenere)
dof_c_tot=dof_c; dof_uv_tot=dof_v; dof_tot=dof;% prima dell'intersezione con i wetnodes
wetnodes = find_wetnodes(elements,cn,g,wdtol,'pred');
drynodes = setdiff(dof_uv_tot,wetnodes);
wetdof = [wetnodes,wetnodes+nov,wetnodes+2*nov]';%veri wetdof
drydof = setdiff(dof_tot,wetdof);%veri drydof
%15/11/2012% dof = wetdof;%veri dof %%% non mi piace, ma almeno non riscrivo tutto...
%15/11/2012% dof_uv_in = intersect(dof_uv_tot,wetnodes); dof_c = intersect(dof_c_tot,wetnodes);
%15/11/2012% drynodes_uv = intersect(dof_uv_tot,drynodes); drynodes_c = intersect(dof_c_tot,drynodes);
ndof=length(dof);ndof_uv_in=length(dof_uv_in);ndof_c=length(dof_c);
toc

disp('risoluzione sistema lineare - Predictor')
% Linear system
    isDir = ~isempty([uDir_in,vDir_in,cDir]);
    isDry = ~isempty(drydof);
    if  isDir | isDry
        if isDir %gestione dei nodi di Dirichlet
            rhs = rhs(dof,1) - aglo(dof,Dnodes)*[uDir_in,vDir_in,cDir,cDir_in]';
            una(nodes0,1) = uDir_in(1:length(nodes0));
            una(nodes1,1) = uDir_in(length(nodes0)+1:length(DirDof_for_uv_in));
            vna(nodes0,1) = vDir_in(1:length(nodes0));
            vna(nodes1,1) = vDir_in(length(nodes0)+1:length(DirDof_for_uv_in));
%            vna(DirDof_for_c,1) = vDir;
            cna(DirDof_for_c,1) = cDir(1:length(DirDof_for_c));
%            cna(nodes0,1) = cDir_in;
        end
        if isDry %gestione dei nodi dry: poniamo tutto a 0
            rhs(drydof) = 0;
            %una(drynodes_uv,1) = 0;
            %vna(drynodes_uv,1) = 0;
            %cna(drynodes_c,1) = 0;
% disp('matrice aglo');
% aglo
figure(102), spy(aglo); title('matrice con Dirichlet')
            aglo(drydof,:) = identita(drydof,:);
            aglo(:,drydof) = identita(:,drydof);
%             aglo(drydof,drydof) = 1.e+10*aglo(drydof,drydof);
% disp('matrice con imposizione su drydof');
% aglo
figure(103); spy(aglo); title('matrice con imposizione su drydof');
% disp('matrice che risolveremo');
% aglo(dof,dof)
figure(104); spy(aglo(dof,dof)); title('matrice che risolveremo');

        end
        temp = aglo(dof,dof)\rhs; % system solution
        una(dof_uv_in,1) = temp(1:ndof_uv_in,1);
        vna(dof_uv_in,1) = temp(ndof_uv_in+1:2*ndof_uv_in,1);
        cna(dof_c,1) = temp(2*ndof_uv_in+1:end,1);
    else
        temp = aglo\rhs;
        una(dof_uv_in,1) = temp(1:ndof_uv_in,1);
        vna(dof_uv_in,1) = temp(ndof_uv_in+1:2*ndof_uv_in,1);
        cna(dof_c,1) = temp(2*ndof_uv_in+1:end,1);
    end

disp('max(u) max(v) max(c)')
[max(una),max(vna),max(cna)]
disp('min(u) min(v) min(c)')
[min(una),min(vna),min(cna)]
pause;

% Aggiornamento variabili altezza e plot Predictor
wdna = cna.^2/4./g;
hna = wdna + h0;
figure(1000);
% pdeplot(vertices,boundaries,elements,'xydata',hna,'contour','on'), axis equal
pdesurf(vertices,elements,hna)
title(strcat('t = ',num2str(t), ' - Predictor'))
pause;
%@<%

% Corrector

[aglo,rhs] = assem_mat_vect_gio_stab_adapt(vertices,elements,boundaries,g,dt,un,vn,cn,una,vna,cna,theta,sigma_res,h0,h_k,f1,f2,f3,wdn,1,t...
    ,unoold,vnoold,cnoold,nx_nodes,ny_nodes); 

%if ~isempty(nodesxy)
%     aglo(nodesxy,:) = sparse(1:nnzxy,nodesxy,ones(1,nnzxy),nnzxy,3*nov,...
%                                  nnzxy);
%     aglo(nodesxy+nov,:) = sparse(1:nnzxy,nodesxy+nov,ones(1,nnzxy),...
%                                      nnzxy,3*nov,nnzxy);
%     rhs([nodesxy,nodesxy+nov]) = 0;   
%end

%@>%

disp('trova nodi wet')
tic
% Find wetnodes (da tenere)
clear wetnodes
wetnodes = find_wetnodes(elements,cn,g,wdtol,'corr');
drynodes = setdiff(dof_uv_tot,wetnodes);
wetdof = [wetnodes,wetnodes+nov,wetnodes+2*nov]';%veri wetdof
drydof = setdiff(dof_tot,wetdof);%veri drydof
%15/11/2012% dof = wetdof;%veri dof %%% non mi piace, ma almeno non riscrivo tutto...
%15/11/2012% dof_uv_in = intersect(dof_uv_tot,wetnodes); dof_c = intersect(dof_c_tot,wetnodes);
%15/11/2012% drynodes_uv = intersect(dof_uv_tot,drynodes); drynodes_c = intersect(dof_c_tot,drynodes);
ndof=length(dof);ndof_uv_in=length(dof_uv_in);ndof_c=length(dof_c);
toc

disp('risoluzione sistema lineare - Corrector')
% Linear system
    isDir = ~isempty([uDir_in,vDir_in,cDir]);
    isDry = ~isempty(drydof);
    if  isDir | isDry
        if isDir %gestione dei nodi di Dirichlet
            rhs = rhs(dof,1) - aglo(dof,Dnodes)*[uDir_in,vDir_in,cDir,cDir_in]';
            un(nodes0,1) = uDir_in(1:length(nodes0));
            un(nodes1,1) = uDir_in(length(nodes0)+1:length(DirDof_for_uv_in));
            vn(nodes0,1) = vDir_in(1:length(nodes0));
            vn(nodes1,1) = vDir_in(length(nodes0)+1:length(DirDof_for_uv_in));
%            vn(DirDof_for_c,1) = vDir;
            cn(DirDof_for_c,1) = cDir(1:length(DirDof_for_c));
%            cn(nodes0,1) = cDir_in;
        end
        if isDry %gestione dei nodi dry: poniamo tutto a 0
            rhs(drydof) = 0;
            %una(drynodes_uv,1) = 0;
            %vna(drynodes_uv,1) = 0;
            %cna(drynodes_c,1) = 0;
            aglo(drydof,:) = identita(drydof,:);
            aglo(:,drydof) = identita(:,drydof);
        end
        temp = aglo(dof,dof)\rhs; % system solution
        un(dof_uv_in,1) = temp(1:ndof_uv_in,1);
        vn(dof_uv_in,1) = temp(ndof_uv_in+1:2*ndof_uv_in,1);
        cn(dof_c,1) = temp(2*ndof_uv_in+1:end,1);
    else
        temp = aglo\rhs;
        un(dof_uv_in,1) = temp(1:ndof_uv_in,1);
        vn(dof_uv_in,1) = temp(ndof_uv_in+1:2*ndof_uv_in,1);
        cn(dof_c,1) = temp(2*ndof_uv_in+1:end,1);
    end
    
disp('max(u) max(v) max(c)')
[max(un),max(vn),max(cn)]
disp('min(u) min(v) min(c)')
[min(un),min(vn),min(cn)]
%@<%
     
% Aggiornamento variabili altezza

wdn = cn.^2/4./g;
hn = wdn + h0;

% Salvataggio
save(strcat(save_path,[num2str(tstart+(t-tspan(1))/dt)]),'un','vn','cn');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST PROCESSING
    
    % CALCOLO PORTATE e VOLUME per le verifiche sulla conservazione della
    % massa
    

[ar] = pdetrg(vertices,elements);     
 mod_v=(un.^2+vn.^2).^0.5;
 
% calcolo della portata attraverso la formula di quadratura di Cavalieri-Simpson 
if sum(bc)~=length(bc)
if in>0    
    

edges_in = find(boundaries(3,:) == in); 
wd_ein (1,:) = wdn(boundaries(1,edges_in)) ;
wd_ein (2,:) = wdn(boundaries(2,edges_in)) ;
wd_e = mean(wd_ein);
un_ein (1,:) = un(boundaries(1,edges_in)) ;
un_ein (2,:) = un(boundaries(2,edges_in)) ;
un_e = mean(un_ein);
dist_in = abs(vertices(2,boundaries(2,edges_in))-vertices(2,boundaries(1,edges_in)));
Q_CSin = dist_in/6.*(wd_ein(1,:).*un_ein(1,:) + 4*wd_e.*un_e + wd_ein(2,:).*un_ein(2,:));
Q_in(p1) = sum(Q_CSin);
disp('Q_in')
Q_in(p1)
end
if in1>0
edges_in1 = find(boundaries(3,:) == in1); 
wd_ein1 (1,:) = wdn(boundaries(1,edges_in1)) ;
wd_ein1 (2,:) = wdn(boundaries(2,edges_in1)) ;
wd_e1 = mean(wd_ein1);
vn_ein1 (1,:) = vn(boundaries(1,edges_in1)) ;
vn_ein1 (2,:) = vn(boundaries(2,edges_in1)) ;
vn_e1 = mean(vn_ein1);
dist_in1 = abs(vertices(1,boundaries(2,edges_in1))-vertices(1,boundaries(1,edges_in1)));
Q_CSin1 = dist_in1/6.*(wd_ein1(1,:).*vn_ein1(1,:) + 4*wd_e1.*vn_e1 + wd_ein1(1,:).*vn_ein1(2,:));
Q_in1(p1) = sum(Q_CSin1);
Q_in1(p1);
end
if out>0
edges_out = find(boundaries(3,:) == out);
wd_eout (1,:) = wdn(boundaries(1,edges_out)) ;
wd_eout (2,:) = wdn(boundaries(2,edges_out)); 
wd_eo = mean(wd_eout);
un_eout (1,:) = un(boundaries(1,edges_out)) ;
un_eout (2,:) = un(boundaries(2,edges_out)) ;
un_eo = mean(un_eout);
dist_out = abs(vertices(2,boundaries(2,edges_out))-vertices(2,boundaries(1,edges_out)));

Q_CSout = dist_out/6.*(wd_eout(1,:).*un_eout(1,:) + 4*wd_eo.*un_eo + wd_eout(2,:).*un_eout(2,:));

Q_out(p1) = sum(Q_CSout);
disp('Q_out')
Q_out(p1)
end
% qx_out_l= wdn(nodesout).*un(nodesout);
% qx_out=0.5*(qx_out_l(1:end-1)+qx_out_l(2:end)).*dist_out';



end


Vol_el = (2*wq*(phiq(1,:)'*wdn(it1)' + phiq(2,:)'*wdn(it2)' + phiq(3,:)'*wdn(it3)')).*ar;

Vol(p1)=sum(Vol_el);

% [err(i)] = bilmassa (0.5*(Q_in(i1)+Q_in(i)),0.5*(Q_out(i1)+Q_out(i)),[Vol(i),Vol(i1)],dt);

 Fr = mod_v./(9.81*wdn).^0.5;
% plot(dt:dt:t,err), pause(1)

% Vol_p=(Vol-Vol(1))/Vol(1);% perdita di volume relativa all'istante iniziale

%

% PLOT variabili principali
% 
% pdesurf(vertices,elements,Fr),pause(1)
% pdesurf(vertices,elements,un),pause(1)
% % % % pdesurf(vertices,elements,cn),pause(1)
% pdesurf(vertices,elements,h0)    
% hold on
% [xx1,yy1] = meshgrid(0:1:100,50);
% h00 = -0.01*xx1 +1;
% hh1 = griddata(vertices(1,:),vertices(2,:),cn.^2/4/9.81,xx1,yy1); 
% uu1 = griddata(vertices(1,:),vertices(2,:),un,xx1,yy1); 

% 
% figure(1)
% plot(xx1,h00 + hh1)
% hold on
% plot(xx1,h00)
% pause(0.2)
% hold off
% figure(2)
% plot(xx1,hh1.*uu1*100)
% pause(0.2)
% % ind = max(vn)

% ERRORE L2 per STOKER 

% h = stoker (vertices(1,:),0,t,g,3,1,1.82,2.397,5.14);
% 
% [errL2] = err_L2(vertices,elements,hn,h)
% err(i) = errL2;
% plot(0:dt:t,Vol_p)
% figure(1)
% pdeplot(vertices,boundaries,elements,'xydata',hn,'contour','on'),axis equal
% pause(0.2)

% [hex] = stoker ([-5:0.01:5],0,t,9.81,0.4,0.016,0.114,1.842,2.14);
% [hex] = stoker ([-5:0.01:5],0,t,9.81,3,1,1.848,2.333,5.082)
% % [hex] = stoker ([-5:0.01:5],0,t,9.81,0.4,0.016,0.114,1.842,2.14);
% % [hex] = stoker ([-5:0.01:5],0,t,9.81,0.4,0.1,0.221,1.019,1.864);
% [xx1,yy1] = meshgrid(-5:0.01:5,0.5);
%     wdn = cn.^2/4/9.81;
%     hh1 = griddata(vertices(1,:),vertices(2,:),wdn,xx1,yy1);
%     plot(xx1,hex,xx1,hh1)
%     pause(0.5)
% CONTROLLO DI CONVERGENZA PER PROBLEMI STEADY

% unoold = unold;
% vnoold = vnold;
% cnoold = cnold;
% 
% unold = un;
% vnold = vn;
% cnold = cn;
% 
% 
% u_qq_old = u_qq;
% v_qq_old = v_qq;
% h_qq_old = h_qq;
% 
% u_qq = phiq(1,:)'*un(it1)' + phiq(2,:)'*un(it2)' + phiq(3,:)'*un(it3)';
% v_qq = phiq(1,:)'*vn(it1)' + phiq(2,:)'*vn(it2)' + phiq(3,:)'*vn(it3)';
% h_qq = phiq(1,:)'*hn(it1)' + phiq(2,:)'*hn(it2)' + phiq(3,:)'*hn(it3)';
% 
% ch_tot  = sum(2*wq*((u_qq-u_qq_old).^2+(v_qq-v_qq_old).^2+(h_qq-h_qq_old).^2).*ar);
% sol_tot = sum(2*wq*((u_qq).^2+(v_qq).^2+(h_qq).^2).*ar);
% 
% solution_change(p) = sqrt(ch_tot/sol_tot);

% disp('solution change')
% figure(2)
% solution_change(p)
% if t/dt>1
%     p0=int64(t/dt-1);
% solution_change(p0)-solution_change(p)
% end

% semilogy(dt:dt:t,solution_change), pause(1)


% if t/dt>1
% if solution_change(p)<5*10^-5
%     break
% % elseif abs(solution_change(p0)-solution_change(p)) < eps
% %     break
% end
% end

figure(1000);
pdeplot(vertices,boundaries,elements,'xydata',hn,'contour','on'), axis equal
title(strcat('t = ',num2str(t), ' - Corrector'))
figure(1001);
pdeplot(vertices,boundaries,elements,'xydata',wdn,'contour','on'), axis equal
title(strcat('t = ',num2str(t), ' - Corrector'))
% % attento = max(abs((vn)))
pause(0.2)
% if mod(t/dt,10)==0
%     pause
% end

%     title(['t = ' num2str(t)])
%     pause(0.1)
% h_sx(i)=hn(125);
% figure(2)
% plot(dt:dt:t,h_sx), pause (1)
% pause
% comptime(p)=toc
% save 'Q_out.mat' Q_out

end % end of time loop

% for i=1:t/dt
%     percVol(i)=(Vol(1)-Vol(i))/Vol(1)*100;
% end

% save 'err.mat' err
save(strcat(save_path,'Vol.mat'), 'Vol');
% save 'solchange.mat' solution_change
% save 'Q_in.mat' Q_in

return
