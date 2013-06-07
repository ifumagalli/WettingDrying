function [un,vn,cin,wdin]=our_iniswe2D_gio(x,y,p,t,time,g,h0,save_path)
%INISWE2D Initial conditions for the shallow water problem
%   [U,V,CIN,ETA]=INISWE2D(X,Y,T) computes the velocity field components U and
%   V and the elevation ETA (with respect to an horizontal reference level)
%   at time T in the nodes (X,Y). CIN = 2*sqrt(g*ETA)
%

% Modificato da Bulgarello,Fumagalli

%   Gaussian-hill test case: eta = 0.2*exp(-0.01*((x-50)^2+(y-50)^2))

%CONDIZIONI INIZIALI PER UN VN

un =  +0.*x;
% un = -1/20.*x;%.*(x<0);
% un = ones(size(x));
% for i=1:size(y)
%    if y(i)< -10
%     un(i) = 0;
%    elseif y(i)> 10
%         un(i) = 0;
%    end 
% end

vn = 0 + 0.*x;

% load 120;
% un=un;
% vn=vn;
% hn=hn;
% wdin = hn - h0;

%CONDIZIONI INIZIALI PER WD

% Onda sinusoidale

% wdin=0 + 0.5*cos(2*pi*(x+20)./80) - h0;

% acqua ferma
% hin = 3*ones(size(x));
% wdin =  max(0,hin - h0); 

% Moto uniforme

% wdin = 0.5 + 0.*x;

% Stoker (all wet) / Ritter (partially dry)
nt =  size(t,2);
[xq,yq,wq,phiq] = basis_on_quad2D('P1',10); 
nq = length(wq);
i1 = t(1,:);
i2 = t(2,:);
i3 = t(3,:);
x1 = p(1,i1);
x2 = p(1,i2);
x3 = p(1,i3);
xqq = phiq(1,:)'*x1 + phiq(2,:)'*x2 + phiq(3,:)'*x3;
% wdqq = 1.e-2*ones(nq,nt);     % Stoker
wdqq = zeros(nq,nt);    % Ritter
[i,j] = find(xqq<=0);
wdqq(i,j) = 6;

wd_el = 2*wq*wdqq;
% wdin = pdeprtni(p,t,wd_el);
% hin = wdin+h0;
% hin = pdeprtni(p,t,wd_el);
% wdin = (hin-h0).*(hin>h0);
wdin = pdeprtni(p,t,wd_el);
wdin = wdin.*(wdin>0);
hin = wdin + h0;

% % Like Stoker/Ritter, but with linear/parabolic transition
% hl = 10; hr = 0; xl = -0.01; xr = 0;
% % hin = hl.*(x<xl) + hr.*(x>xr) + ((hl-hr)/(xr-xl).*(xr-x)+hr).*(x>=xl).*(x<=xr);
% % hin = hl.*(x<0) + hr.*(x>xr) + (hl/xr^2.*(x-xr).^2+hr).*(x>=0).*(x<=xr);
% hin = hl.*(x<xl) + hr.*(x>xr) + (x.^2/20).*(x>=xl).*(x<=xr);
% wdin = (hin-h0).*(hin>h0);

%Parabolic bowl
% versione vecchia: non sappiamo cosa fa: non ci torna parabolic...
% X = 1;
% Y = -0.41884;
% alpha = 1.61e-7;
% 
% [h,eta,u,v,z,tau] = parabolic (x,y,0,alpha,X,Y);
% wdin = h;
% 
% pdesurf(p,t,h)
% pause
% hin = 4;
% wdin = 4-h0;
% 
% load 99
% un=un;
% vn=vn;
% hin = hn;
% wdin = hin;

% versione nuova
% hin = 1.e-7*(8-.3*(x.^2+y.^2));
% hin = 0.55-0.0375*(x.^2+y.^2);
% underground_idxs = find(hin<=h0);
% wdin = hin-h0;
% wdin(underground_idxs) = 0;

% PLANE ON PARABOLIC BOWL
% hin = 0.475-0.03*x;
% figure(1000); pdesurf(p,t,h0); hold on; pdesurf(p,t,hin);
% underground_idxs = find(hin<=h0);
% wdin = hin-h0;
% wdin(underground_idxs) = 0;
% title('Initial condition')
% pause

% Radial Dam Break
% 
% c = [0,0];
% nt =  size(t,2);
% % r = ((x-c(1)).^2+(y-c(2)).^2).^0.5;
% [xq,yq,wq,phiq] = basis_on_quad2D('P1',10); 
% nq = length(wq);
% i1 = t(1,:);
% i2 = t(2,:);
% i3 = t(3,:);
% x1 = p(1,i1);
% x2 = p(1,i2);
% x3 = p(1,i3);
% y1 = p(2,i1);
% y2 = p(2,i2);
% y3 = p(2,i3);
% xqq = phiq(1,:)'*x1 + phiq(2,:)'*x2 + phiq(3,:)'*x3;
% yqq = phiq(1,:)'*y1 + phiq(2,:)'*y2 + phiq(3,:)'*y3;
% rqq = sqrt(xqq.^2 + yqq.^2);
% 
% %@>%
% wdqq = zeros(nq,nt);
% [i,j] = find(rqq<.5);
% wdqq(i,j) = 1;
% %@<%
% 
% wd_el = 2*wq*wdqq;
% wdin = pdeprtni(p,t,wd_el);
% % wdin = wdin';
% % for i=1:size(r)
% %    if r(i)>0.5
% %        wdin(i)=1;
% %    else   
% %       wdin(i) = 2;
% %    end
% % end
% % wdin=wdin';

% Onda di marea
% un =  -0.3 +0.*x;
% vn = 0 + 0.*x;
% wdin= 0.4 +0.*x;

%panettone S&H

% wdin=(1+0.2*sin(0.02*pi*(x-25))).*(1+0.2*sin(0.02*pi*(y-25)))
% un=0 + 0.*x;
% vn=0 + 0.*y;
% Metto a 0 eventuali wd<0

% % canale TFDC
% 
% vn = 0 + 0.*x;
% wdin= 1.458-(1.458-0.76).*x/500;
% un = 0 + 0.*x;

% canale allargamento steady

% vn = 0 + 0.*x;
% wdin= 0.5 +0.*x;

% ch1 = find(y<=40);

% OBLIQUE

% wdin = 1 + 0.*x;
% un=8.57 + 0.*x; 

% CONFLUENZA

% un =8 + 0.*x; % conf
% wdin = 0.4 + 0.*x;
% un =2 + 0.*x; % conf lenta
% wdin =1.869 +0.*x % conf lenta


% ponte
% 
% vn = 0 + 0.*x;
% hin = 2.658 - 0.001.*x;
% wdin= hin-h0;
% un =  0 +0.*x;
% % 
% % ponte adapt
% load ini_ponte.mat
% 
% un = un;
% vn = vn;
% hin = hn;
% 
% wdin = hin - h0;

% 
% neg = find(wdin<0);
% wdin(neg)=0;
%  xin =1 + 5*exp(-0.25*((x-25).^2+(y-25).^2));% Gaussiana simona
%  xin= +0.5+0.001.*x; %altezza costante
%   xin =0.5+exp(-0.001*((x).^2+(y).^2))

% Cilindro
% un = sqrt(11.76)*0.2*(sech(sqrt(0.15)*(x+20))).^2./(5+0.2*(sech(sqrt(0.15)*(x+20))).^2);
% vn = 0 + 0*x;
% xin = +0.2*(sech(sqrt(0.15)*(x+20))).^2;

% Nozzle
% un = sqrt(11.76)*0.2*(sech(sqrt(0.15)*(x-7))).^2./(5+0.2*(sech(sqrt(0.15)*(x-7))).^2);

% xin = 5+0.2*(sech(sqrt(0.15)*(x-7))).^2;

% sudden

% wdin = 2 +0.*x -h0;



% BUMP

% hin = 2 +0.*x;
% wdin= hin-h0; 
% un = 0 +0.*x;
% vn = 0 +0.*x;

% HYPPIE

% hin = 1 - 0.001*x;
% wdin= hin-h0; 
% un = x.^2 - x +0.*x;
% vn = y.^2 - y +0.*x;

% Uniforme

% un = 1.17 + 0.*x;
% wdin = 0.05 + (0.6-0.05)/100.*x;
% wdin = 0.38 + 0.*x

% pv = find (x<50);
% pl = find (x>=50);
% wdin(pv) = 0.38;
% wdin(pl) = 1;
% wdin = wdin';

% 
% i011 = find(x>=0);
% i04 = find(x<0);
% hin(i04) = 0.4;
% hin(i011) = 0.016;
% % hin = hin';
% % hin = max(hin,h0);
% wdin = hin' - h0;
% cin = 2*(g*wdin).^0.5;

cin = 2*(g*wdin).^0.5;
% % 
% cin = 10 + 2*cos(x)+ tanh(10*(x-y));
% cin = 10 + y.^2 + x.^2;
% cin = 10 + 0.*x;
% wdin = cin.^2./(4*g);
cn=cin;
save(strcat(save_path,'init.mat'),'un','vn','cn');
pfig=figure(1000); pdesurf(p,t,h0); hold on; pdesurf(p,t,hin);
title('Initial condition')
print(pfig,'-deps',strcat(save_path,'AAAinitial_both','.eps'));
print(pfig,'-djpeg',strcat(save_path,'AAAinitial_both','.jpeg'));
saveas(pfig,strcat(save_path,'AAAinitial_both','.fig'));
pause

return
