function [u_ex,h_ex,shockvelocity_ex] = exact_Stoker_Ritter(vertices,elements,boundaries,save_path,tspan,h0,hl,hr,g,n,do_print)

%   [u_ex,h_ex,shockvelocity_ex] = exact_Stoker_Ritter(vertices,elements,boundaries,save_path,tspan,h0,hl,hr,g)
% 
% Computes exact solution of Stoker/Ritter problem with initial heigth
% discontinuity on x=0 and initial velocity u=0.
%
% tspan:        time interval in which to compute the solution
% h0:           bathimetry
% hl, hr:       initial total elevation for x<0, x>0, respectively
% g:            gravity acceleration
% n:            Manning friction coefficient
% u_ex,h_ex:    exact solution: there is 0 wherever there's no water
% do_print:     =1 if you want to have some plots (default =0)
% NB: v_ex = 0
% 
% [solution from a combination of:
% O. Delestre, C. Lucas, P.-A. Ksinant, F. Darboux, C. Laguerre, T. N. T. Vo, F. James, and S. Cordier, “SWASHES : a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies.” pp. 1–40, 2012.
% J. Hudson, “Numerical Techniques for Morphodynamic Modelling,” Reading, 2001.

N = max(size(vertices));
t=linspace(tspan(1),tspan(2),tspan(3)+1);

if(nargin == 9)
    do_print = 0;
end

isRitter = (hr==0);

shockvelocity_ex = [];
if(n==0)
    syms cm
    % cm = ~isRitter * solve(-8*g*hr*cm^2*(g*hl-cm^2)^2+(cm^2-g*hr)^2*(cm^2+g*hr),cm);
    cm = ~isRitter * roots([1, 0, -9*g*hr, 16*g*hr*sqrt(g*hl), -g^2*hr*(hr+8*hl), 0, (g*hr)^3]);
        % NB: l'equazione che risolviamo per trovare cm non è quella di questo
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

    if isRitter
        shockvelocity_ex = ex_xB(1);
    else
        shockvelocity_ex = ex_xC(1);
    end

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
    
else
    C = 1/n;%*hl^(1/6);
    ex_xA = @(t) -t*sqrt(g*hl);
    ex_xB = @(t) t*2*sqrt(g*hl);
    
    alpha1 = @(x,t) 6/5*(1./(2-x*(1./(t.*sqrt(g*hl)))))-2/3+4/135*sqrt(3)*(2-x*(1./(t.*sqrt(g*hl)))).^(3/2);
    alpha2 = @(x,t) 12*(1./(2-x*(1./(t.*sqrt(g*hl)))))-8/3+8/189*sqrt(3)*(2-x*(1./(t.*sqrt(g*hl)))).^(3/2)-108/7*(1./(2-x*(1./(t.*sqrt(g*hl))))).^2;
    ex_hco = @(x,t) 1/g*(2/3*sqrt(g*hl)-x*(1./(3*t))+g^2*(1./C.^2).*alpha1(x,t).*t).^2;
    ex_uco = @(x,t) 2/3*(sqrt(g*hl))+2*x*(1./(3*t))+g^2*(1./C.^2).*alpha2(x,t).*t;

    [utip_ex,idx] = max(ex_uco(vertices(1,:),t).*(vertices(1,:)<ex_xB(t)).*(vertices(1,:)>ex_xA(t)));
    xT_ex = vertices(1,idx);
    
%     ex_h = @(x,t) hl*(x*ones(size(t))<=ex_xA(ones(size(x))*t)) ...
%             + ex_hco(x,t).*(x*ones(size(t))<ex_xB(ones(size(x))*t)).*(x*ones(size(t))>ex_xA(ones(size(x))*t)) ...
%             + hr*(x*ones(size(t))>=ex_xB(ones(size(x))*t));
        
    ex_h = @(x,t) hl*(x*ones(size(t))<=ex_xA(ones(size(x))*t)) ...
            + ex_hco(x,t).*((x<xT_ex)*ones(size(t))).*(x*ones(size(t))>ex_xA(ones(size(x))*t)) ...
            + hr*(x>=xT_ex);
        
    ex_u = @(x,t) 0*(x*ones(size(t))<=ex_xA(ones(size(x))*t)) ...
            + ex_uco(x,t).*((x<xT_ex)*ones(size(t))).*(x*ones(size(t))>ex_xA(ones(size(x))*t)) ...
            + utip_ex*(x*ones(size(t))<ex_xB(ones(size(x))*t)).*((x>=xT_ex)*ones(size(t))) ...
            + 0*(x*ones(size(t))>=ex_xB(ones(size(x))*t));

    h_ex = ex_h(vertices(1,:)',t);
    u_ex = ex_u(vertices(1,:)',t);
    shockvelocity_ex = ex_xB(1);
end

if do_print==1

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
    idxs_ymin = find(vertices(2,:) == ymin);
    ymax = max(vertices(2,:));
    idxs_ymax = find(vertices(2,:) == ymax);
    [~,idx] = min(abs(vertices(2,:)-(ymax+ymin)/2));
    ymed = vertices(2,idx);
    idxs_ymed = find(vertices(2,:) == ymed);

    for i=1:length(t)
    
        pfig = figure(501);  set(gcf,'Visible','off');
        pdeplot(vertices,boundaries,elements,'zdata',h_ex(:,i),'zstyle','discontinuous')
        title(strcat('exact h (where wet) at t=',num2str(t(i))))%,'  Vol=',num2str(volume)));
        % print(pfig,'-deps',strcat(save_path,'h',num2str(t(i),'%.3f'),'exact','.eps'));
        print(pfig,'-djpeg',strcat(save_path,'h',num2str(t(i),'%.3f'),'exact','.jpeg'));
        saveas(pfig,strcat(save_path,'h',num2str(t(i),'%.3f'),'exact','.fig'));

        pfig = figure(502);  set(gcf,'Visible','off');
        plot(vertices(1,idxs_ymin),h_ex(idxs_ymin,i),vertices(1,idxs_ymax),h_ex(idxs_ymax,i),vertices(1,idxs_ymed),h_ex(idxs_ymed,i))
        legend('y = y_m_i_n','y = y_m_a_x','y = y_m_e_d')
        title(strcat('h_e_x at t = ',num2str(t(i))));
        % print(pfig,'-deps',strcat(save_path,'hn_wet',num2str(t,'%.3f'),'apred','.eps'));
        print(pfig,'-djpeg',strcat(save_path,'h_ysection',num2str(t(i),'%.3f'),'exact','.jpeg'));
        saveas(pfig,strcat(save_path,'h_ysection',num2str(t(i),'%.3f'),'exact','.fig'));

        % close all
        pause(0.1);
    end
end

% max(max(h_ex))

end