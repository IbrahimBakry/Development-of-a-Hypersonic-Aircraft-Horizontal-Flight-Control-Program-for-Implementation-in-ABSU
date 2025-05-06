% ********************************************************
% Function Aerodynamics
% Provision of Aerodynamic Forces and Moments in Body Fixed 
% Coordinates
% ********************************************************

function [F,F_v,M,M_v,MG1,rho,g,ma] = aerodynamics(h,Vabs,rg1,Mg,alfa,beta)

[m,Sref,Lref,bref,rhos,gs,bet,omegae,rearth] = constants;

arc = pi/180.;

% ********************************************************
% Generation of Aerodynamics of Space Vehicle
% ********************************************************

% rho = rhos*exp(-bet*h);
H = 10351.8-(0.0368512).*h-(1.02368e-5).*h.^2+(2.63363e-10).*h.^3;
rho = rhos.*exp(-h./H);
g   = gs*(rearth/rg1)^2;

norm = rho/2.*Vabs^2*Sref;
norm2 = rho.*Vabs*Sref;

% algr = alfa/arc;
% begr = beta/arc;

% ********************************************************
% X-38 Longitudinal and Lateral Aerodynamics
% ********************************************************

CL =@(alfa) -0.0057 + 0.8671*alfa + 3.5556*alfa^2;
CD =@(alfa) -0.0101 + 0.1031*alfa + 2.1353*alfa^2;
Cm =@(alfa)  0.1116 + 0.1384*alfa - 0.7629*alfa^2;
CQ =@(beta)  0.0000 - 0.1038*beta + 0.0299*beta^2;
Cn =@(beta)  0.0000 + 0.1080*beta - 0.0120*beta^2;
Cll = 0.;

D = norm*CD(alfa);
Y = norm*CQ(beta);
L = norm*CL(alfa);
F = [D Y L];
F_v = [norm2*CD(alfa) norm2*CQ(beta) norm2*CL(alfa)];

ma(1) = Cll*Lref;
ma(2) = Cm(alfa)*Lref;
ma(3) = Cn(beta)*Lref;

M = norm*ma';
M_v = [norm2*ma(1) norm2*ma(2) norm2*ma(3)];

G(1) = 0.;
G(2) = 0.;
G(3) = m*g;

MG1=Mg*G';

% ********************************************************
% End of Function Aerodynamics 
% ********************************************************
