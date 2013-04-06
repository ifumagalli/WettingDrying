function [u_ex,h_ex] = exact_Stoker_Ritter(vertices,elements,boundaries,save_path,tspan,h0,hl,hr,g,uin,vin)

%       [u_ex,h_ex] = exact_Stoker_Ritter(vertices,elements,boundaries,save_path,tspan,h0,hl,hr,g,uin,vin)
% 
% Computes exact solution of Stoker/Ritter problem with initial
% discontinuity on x=0.
%
% h0:           bathimetry
% hl, hr:       initial total elevation for x<0, x>0, respectively
% g:            gravity acceleration
% uin, vin:     initial x and y velocity distributions (default: uin=vin=0)
% u_ex,h_ex: exact solution: there is 0 wherever there's no water
% NB: v_ex = 0
% 
% [solution from: Delestre, O., Lucas, C., Ksinant, P.-A., Darboux, F.,
% Laguerre, C., Vo, T. N. T., James, F., et al. (2012). SWASHES : a
% compilation of Shallow Water Analytic Solutions for Hydraulic and
% Environmental Studies.]

N = max(size(vertices));
t=linspace(tspan(1),tspan(2),tspan(3)+1);
if(nargin == 9)
    uin = zeros(N,1);
    vin = uin;
end

isRitter = (hr==0);

syms cm
% cm = ~isRitter * solve(-8*g*hr*cm^2*(g*hl-cm^2)^2+(cm^2-g*hr)^2*(cm^2+g*hr),cm);
cm = ~isRitter * roots([1, 0, -9*g*hr, 16*g*hr*sqrt(g*hl), -g^2*hr*(hr+8*hl), 0, (g*hr)^3])
    % NB: l'equazione che risolviamo per trovare cm non � quella di questo
    % articolo succitato, ma l'abbiamo ricavata da un misto di questo
    % articolo e l'Hudson, ponendo cm=c2
    % PROVARE CON L'EQUAZIONE DELLO SWASHES
cm = max(cm(find(imag(cm)==0)))
ex_hraref = @(x,t) 4/(9*g)*(sqrt(g*hl)-x*(1./(2*t))).^2;
ex_hm = @(x,t) ~isRitter * cm^2/g*ones(size(x))*ones(size(t));
ex_uraref = @(x,t) 2/3*(sqrt(g*hl)+x*(1./t));
ex_um = @(x,t) ~isRitter * 2*(sqrt(g*hl)-cm);

ex_xA = @(t) -t*sqrt(g*hl);
ex_xB = @(t) t*(2*sqrt(g*hl)-3*cm);
ex_xC = @(t) t*2*cm^2*(sqrt(g*hl)-cm)/(cm^2-g*hr);

ex_h = @(x,t) hl*(x*ones(size(t))<=ex_xA(ones(size(x))*t)) ...
        + ex_hraref(x,t).*(x*ones(size(t))<ex_xB(ones(size(x))*t)).*(x*ones(size(t))>ex_xA(ones(size(x))*t)) ...
        + ex_hm(x,t).*(x*ones(size(t))<ex_xC(ones(size(x))*t)).*(x*ones(size(t))>=ex_xB(ones(size(x))*t)) ...
        + hr*(x*ones(size(t))>=ex_xC(ones(size(x))*t));
    
ex_u = @(x,t) 0*(x*ones(size(t))<=ex_xA(ones(size(x))*t)) ...
        + ex_uraref(x,t).*(x*ones(size(t))<ex_xB(ones(size(x))*t)).*(x*ones(size(t))>ex_xA(ones(size(x))*t)) ...
        + ex_um(x,t).*(x*ones(size(t))<ex_xC(ones(size(x))*t)).*(x*ones(size(t))>=ex_xB(ones(size(x))*t)) ...
        + 0*(x*ones(size(t))>=ex_xC(ones(size(x))*t));

h_ex = ex_h(vertices(1,:)',t);
u_ex = ex_u(vertices(1,:)',t);
    
% figure(500)
% pdeplot(vertices,boundaries,elements,'xydata',ex_h(vertices(1,:),0),'contour','on')
% figure(501)
% pdesurf(vertices,elements,ex_h(vertices(1,:),0)')
% % hold on
% % pdesurf(vertices,elements,hin(vertices(1,:),0)'-h0(vertices(1,:),0)')
% hold on
% pdesurf(vertices,elements,hin(vertices(1,:),0)')
% % figure(502)
% % pdesurf(vertices,elements,hin(vertices(1,:)',0)-ex_h(vertices(1,:)',0))

% mesh utilities
ymin = min(vertices(2,:));
idxs_ymin = find(vertices(2,:) == ymin)
ymax = max(vertices(2,:));
idxs_ymax = find(vertices(2,:) == ymax);
[~,idx] = min(abs(vertices(2,:)-(ymax+ymin)/2));
ymed = vertices(2,idx);
idxs_ymed = find(vertices(2,:) == ymed);

for i=1:length(t)
    
    pfig = figure(501);
    pdeplot(vertices,boundaries,elements,'zdata',h_ex(:,i),'zstyle','discontinuous')
    title(strcat('exact h (where wet) at t=',num2str(t(i))))%,'  Vol=',num2str(volume)));
    saveas(pfig,strcat(save_path,'h',num2str(t(i),'%.3f'),'exact','.fig'));

    pfig = figure(502);
    plot(vertices(1,idxs_ymin),h_ex(idxs_ymin,i),vertices(1,idxs_ymax),h_ex(idxs_ymax,i),vertices(1,idxs_ymed),h_ex(idxs_ymed,i))
    legend('y = y_m_i_n','y = y_m_a_x','y = y_m_e_d')
    title(strcat('h_e_x at t = ',num2str(t(i))));
    % print(pfig,'-deps',strcat(save_path,'hn_wet',num2str(t,'%.3f'),'apred','.eps'));
    print(pfig,'-djpeg',strcat(save_path,'h_ysection',num2str(t(i),'%.3f'),'exact','.jpeg'));
    saveas(pfig,strcat(save_path,'h_ysection',num2str(t(i),'%.3f'),'exact','.fig'));

    % close all
    pause(0.1);
end

max(max(h_ex))

end