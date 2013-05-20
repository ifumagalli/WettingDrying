function [cdir,cdir1,vdir]=xidirswe2D(x,y,boundaries,vertices,h0,t,g,nlato,nlato1)
%XIDIRSWE2D Dirichlet condition for the elevation in the shallow water problem
%   ETA=XIDIRSWE2D(X,Y,T) assign the elevation ETA 
%   (with respect to an horizontal reference level) at time T 
%   in the nodes (X,Y)
%

% xi1 = 0 + 0.*x;
% nbn = 2;
vdir=[];
% if nlato1 == 0 

% MAREA


% wddir = 0.4 + 0.8*(1+cos(2*pi*t/44700-pi))+0.*x;
% wddir = 0.4 +0.*x;% + 0.8*(1+cos(2*pi*t/44700-pi))+0.*x; 
% RISALTO

% wddir = 1 + 0.*x; 

%OBS

% hdir = 2.4 - x*0.001;
% wddir = hdir-h0'; 

% xi= wddir' + h0;%elevazione

% SE
% wddir = 2 +0.*x; % tolta
% 

% conf lenta
% wddir = 1.869+0.*x;

% conf veloce

% wddir = 0.4 +0.*x;

% % oblique
% wddir = 1 + 0.*x;


% bump 
% subcritico
% hdir = 2 + 0.*x;
% wddir = 2 - h0';
% % 
% % 
% cdir=2*(g*wddir).^0.5 ;
% 
vdir=0+0.*x;
cdir1 =[];
% 
% cdir = 10 + 2*cos(x) + tanh(10*(x-y));

% Moto uniforme
wddir = 1;
% if t<100
%   wddir = 0.38+0.0062*t;
%   else
%       wddir = 1;
% end
cdir = 2*sqrt(g*wddir) + 0.*x;
% % cdir = 10 + x.^2;
% else
% 
% edges0 = find(boundaries(3,:) == nlato); 
% nodes0 = unique(boundaries(1:nbn,edges0));
% 
% wddir(1:length(nodes0)) = 1.869;
% 
% cdir=2*(g*wddir).^0.5;
% edges1 = find(boundaries(3,:) == nlato1); 
% nodes1 = unique(boundaries(1:nbn,edges1));
% 
% 
% wddir1(1:length(nodes1)) = 3 - h0(nodes1)';
% 
% cdir1=2*(g*wddir1).^0.5;
% end

return
