function [uDir_in0,vDir_in0,uDir_in1,vDir_in1,cDir_in] = uvdir_inswe2D(x,y,t,wd,nodes,g,nlato,nlato1,boundaries,vertices)
%UVDIRSWE2D Dirichlet condition for the velocity field
% uvdir_inswe2D(vertices(1,DirDof_for_uv_in),vertices(2,DirDof_for_uv_in),t,wdn,DirDof_for_uv_in,g,in,in1,boundaries,vertices);
tic
if nlato == 0
    vDir_in1 =[];
    uDir_in1 =[];
    vDir_in0 =[];
    uDir_in0 =[];
    cDir_in =[];
else
    % Q=250;
    [nr,nc]=size(y);
    % q_x=Q/(nc-1);
    % q_y=0;

    % uDir_in=(q_x+0.*y)./wd(nodes,1)';

    %Risalto

    % uDir_in = -10 +0.*x;
    % uDir_in = 4.15 + 0.*x;


    % uDir_in = 0.*x;

    % cDir_in = 2*sqrt(0.38*9.81) + 0.*x;
    % risalto
    % wdDir_in = 0.3 + 0.*x; %
    % cDir_in=2*(g*wdDir_in).^0.5; % lascio attiva anche negli altri casi anche se poi non la uso
    % 
    % % canale TFDC
    % 
    % Q = 80;

    % allargamento steady
    %  if t<18
    %     if t < 6 
    % Q = 10*t;
    % end
    % if t >= 6 & t < 12
    % Q = 60;
    % end
    % if t >= 12 & t < 18
    % Q = 60 - (t - 12)*10;
    % end
    % else 
    %     Q = 0;
    % end
    % % % % if t<10
    % % % % Q = 22.1; % bump
    % % % % Q = 60;   % conf
    % % % % Q = 10*t % Sudden
    % % % % Q = 100;
    % % % % else 
    % % % % Q = 257.1; % oblique
    % % % % Q = 100 ; % Sudden    
    % % % % Q = 22.1; % bump
    % % % % Q = 60;   % conf
    % % % % Q = 100;
    % % % % end
    %Marea
    Q = 0; 
    nbn = 2;
    edges = find(boundaries(3,:) == nlato); 
    nodes0 = unique(boundaries(1:nbn,edges));
    % dist = abs(vertices(2,nodes0)-vertices(2,nodes0));
    % b = sum(dist);
    b = max(vertices(2,nodes0))-min(vertices(2,nodes0));


    q = Q/b;
    % q = 0;

    uDir_in0(1:length(nodes0)) = q./(wd(nodes0));
    % uDir_in0(1:length(nodes0))=-0.3;
    % % uDir_in0(1:length(nodes0)) = 8.57; % oblique
    vDir_in0(1:length(nodes0)) = 0;
    cDir_in =[]; 
%     cDir_in(1:length(nodes0)) = 2*sqrt(g*3);
    if nlato1>0
        if t<20
        Q1 = t;
        else
            Q1 = 20;
        end

        nbn = 2;
        edges1 = find(boundaries(3,:) == nlato1);
        nodes1 = unique(boundaries(1:nbn,edges1));
        dist1 = abs(vertices(1,boundaries(2,edges1))-vertices(1,boundaries(1,edges1)));
        b1 = sum(dist1);

        q1 = Q1/b1;

        vDir_in1(1:length(nodes1)) = -q1./wd(nodes1);

        uDir_in1(1:length(nodes1)) = 0;
        else
        vDir_in1 =[];
        uDir_in1 =[];
    end
    % % curva inflow con profilo parabolico
    % uDir_in = 0.5*(-y.^2+b*y)/100;
    % 
    % uDir_in = 1-exp(-(-y.^2+b*y));


    % q = 5/sum(dist);
    % % q_e = q.*dist
    % q_e =q+0.*x;
    % 
    % % u_e = q_e'./(0.5*(wd(nodes(1:end-1),1)+wd(nodes(2:end))));
    % 
    % q_n = 0.5*([q_e(1:end)]+[q_e(1:end)]);
    %
    % u_n = q_n'./wd(nodes);
    % 
    % % uDir_in(1:length(nodes)) = Q/sum(0.5*(wd(nodes(1:end-1),1)+wd(nodes(2:end),1))'.*dist);
    % uDir_in = u_n';



    % sudden
    % 
    % if t < 10
    %     Q = 20*t;
    % else
    %     Q=200;
    % end
    % 
    % uDir_in = Q/sum(wd(nodes(1:end-1),1)'.*dist)+0.*x;
    % vDir_in = 0 + 0.*x;


    % [nr,nc]=size(y);
    % 
    % q_y=0;
     toc
    % uDir_in0 = 0+ 0.*x;
    % vDir_in0 = y.^2 - y;
end
return
