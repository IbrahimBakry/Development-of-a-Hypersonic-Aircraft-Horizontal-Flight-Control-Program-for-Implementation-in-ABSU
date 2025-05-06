function [D,Y,L,Dv,Yv,Lv,g,rho] = aero1(in)

Sref   = 21.672;
rhos   = 1.225;
gs     = 9.80665;
bet    = 0.000140845;
rearth = 6.38E6;


arc = pi/180.

Vabs = in(1)
h = in(2)
alfa = in(3)

% ********************************************************
% Generation of Aerodynamics of Space Vehicle
% ********************************************************

R = h + rearth;
beta = 0;

rho = rhos*exp(-bet*h);
% H = 10351.8-(0.0368512).*h-(1.02368e-5).*h.^2+(2.63363e-10).*h.^3;
% rho = rhos.*exp(-h./H);
g   = gs*(rearth/R)^2;

norm = rho/2.*Vabs^2*Sref;
norm2 = rho.*Vabs*Sref;

% algr = alfa/arc;
% begr = beta/arc;

% ********************************************************
% X-38 Longitudinal and Lateral Aerodynamics
% ********************************************************

CL = -0.0057 + 0.8671*alfa + 3.5556*alfa^2;
CD = -0.0101 + 0.1031*alfa + 2.1353*alfa^2;
CQ =  0.0000 - 0.1038*beta + 0.0299*beta^2;

D = norm*CD;
Y = norm*CQ;
L = norm*CL;
F = [D Y L];
Dv = norm2*CD;
Yv = norm2*CQ;
Lv = norm2*CL;
F_v = [Dv Yv Lv];


% ********************************************************
% End of Function Aerodynamics 
% ********************************************************
