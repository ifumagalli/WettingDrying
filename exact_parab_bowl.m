function [u_ex,v_ex,h_ex] = exact_parab_bowl(vertices,elements,boundaries,tspan,a,b,c,g,uin,vin)

%       exact_parab_bowl(vertices,elements,boundaries,tspan,a,b,c,g)
%       exact_parab_bowl(vertices,elements,boundaries,tspan,a,b,c,g,uin,vin)
% 
% Computes exact solution of free oscillations in a paraboloidic basin
% from initial conditions given by the data that are passed to the function:
% a concave paraboloidic free surface on a convex paraboloidic bowl, both
% with vertex in (x,y) = (0,0).
% NB: The mesh (given through vertices,elements,boundaries) is supposed to
% be a square with the centre in (0,0) and the basin has to have floor at
% z=0 on his vertex (0,0)
% 
% Basin:        h0 = a*(x^2+y^2);
% Free surface: h = h0 + wd = c - b*(x^2+y^2);
% g:            gravity acceleration
% uin, vin:     initial x and y velocity distributions (default: uin=vin=0)
% u_ex,v_ex,h_ex: exact solution: there is 0 wherever there's no water
% 
% [solution, adapted for a different input, from: Andrea Balzano,
% Evaluation of methods for numerical simulation of wetting and drying in
% shallow water flow models, Coastal Engineering 34 (1998) 83-107]

N = max(size(vertices));
t=linspace(tspan(1),tspan(2),tspan(3)+1);
if(nargin == 6)
    uin = zeros(N,1);
    vin = uin;
end

% a,b,c -> r0,L,D0,h0
r0 = sqrt(c/(a+b)); % initial radius of shoreline
% L = ( 3.0*r0/(2.0*a) * (c-r0^2*(a+b)/3.0) )^(1.0/3.0);
% L = ( 3*r0^2*c/(2*a) ) ^0.25
L = (r0^2*c/a)^0.25
    % radius of the basin at the equilibrium level z_eq = a*L^2;
D0 = a*L^2; % initial water elevation in (0,0) w.r.t. equilibrium level
h0 = c-a*L^2; % initial total elevation in (0,0) w.r.t. equilibrium level

omega = sqrt( 8*g*D0/L^2 );
T=2*pi/omega
A = ( (D0 + h0)^2 - D0^2 ) / ( (D0 + h0)^2 + D0^2 );

ex_u = @(x,y,t) 0.5*omega*x.*A * (sin(omega*t) ./ (1 - A*cos(omega*t)));
ex_v = @(x,y,t) 0.5*omega*y.*A * (sin(omega*t) ./ (1 - A*cos(omega*t)));
ex_h = @(x,y,t) a*L^2 + D0*( sqrt(1-A^2)./(ones(size(x))*(1-A*cos(omega*t))) ...
                - 1 - (x.^2+y.^2)/L^2 * ((1-A^2)./(1-A*cos(omega*t)).^2 - 1 ) );
            
hin = @(x,y,t) 0.55-0.0375*(x.^2+y.^2)*ones(size(t));
h0 = @(x,y,t) 0.1*(x.^2+y.^2)*ones(size(t));
disp('eeeee')
max(max(hin(vertices(1,:),vertices(2,:),0)-h0(vertices(1,:),vertices(2,:),0)-ex_h(vertices(1,:),vertices(2,:),0)))
min(min(hin(vertices(1,:),vertices(2,:),0)-h0(vertices(1,:),vertices(2,:),0)-ex_h(vertices(1,:),vertices(2,:),0)))
max(max(hin(vertices(1,:),vertices(2,:),0)-ex_h(vertices(1,:),vertices(2,:),0)))
min(min(hin(vertices(1,:),vertices(2,:),0)-ex_h(vertices(1,:),vertices(2,:),0)))
a*L^2>c

figure(500)
size(ex_h(vertices(1,:),vertices(2,:),0))
pdeplot(vertices,boundaries,elements,'xydata',ex_h(vertices(1,:),vertices(2,:),0),'contour','on')
figure(501)
pdesurf(vertices,elements,ex_h(vertices(1,:),vertices(2,:),0)')
% hold on
% pdesurf(vertices,elements,hin(vertices(1,:),vertices(2,:),0)'-h0(vertices(1,:),vertices(2,:),0)')
hold on
pdesurf(vertices,elements,hin(vertices(1,:),vertices(2,:),0)')
% figure(502)
% pdesurf(vertices,elements,hin(vertices(1,:)',vertices(2,:)',0)-ex_h(vertices(1,:)',vertices(2,:)',0))

A=[ones(max(size(vertices)),1),vertices(1,:)',vertices(2,:)',vertices(1,:)'.^2,vertices(2,:)'.^2];
b=hin(vertices(1,:)',vertices(2,:)',0)-h0(vertices(1,:)',vertices(2,:)',0)-ex_h(vertices(1,:)',vertices(2,:)',0);
coeff=A\b

b=hin(vertices(1,:)',vertices(2,:)',0)-ex_h(vertices(1,:)',vertices(2,:)',0);
coeff=A\b

u_ex = ex_u(vertices(1,:)',vertices(2,:)',t).*(ex_h(vertices(1,:)',vertices(2,:)',t)>h0(vertices(1,:)',vertices(2,:)',t));
v_ex = ex_v(vertices(1,:)',vertices(2,:)',t).*(ex_h(vertices(1,:)',vertices(2,:)',t)>h0(vertices(1,:)',vertices(2,:)',t));
h_ex = ex_h(vertices(1,:)',vertices(2,:)',t).*(ex_h(vertices(1,:)',vertices(2,:)',t)>h0(vertices(1,:)',vertices(2,:)',t));

figure(503)
pdeplot(vertices,boundaries,elements,'zdata',h_ex(:,1),'zstyle','discontinuous')

pause;
close all
figure(503)
for i=1:length(t)
%     volume = h_ex(:,i)-h0(vertices(1,:)',vertices(2,:)',0);
    pdeplot(vertices,boundaries,elements,'zdata',h_ex(:,i),'zstyle','discontinuous')
    zlim([0 0.6])
    title(strcat('exact h (where wet) at t=',num2str(t(i))))%,'  Vol=',num2str(volume)));
    pause(0.1);
end

max(max(h_ex))

end