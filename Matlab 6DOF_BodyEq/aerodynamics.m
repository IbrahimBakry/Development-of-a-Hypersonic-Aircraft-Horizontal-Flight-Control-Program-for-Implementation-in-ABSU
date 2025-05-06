% ********************************************************
% Function Aerodynamics
% Provision of Aerodynamic Forces and Moments in Body Fixed 
% Coordinates
% ********************************************************

function [A,M,MG1,CL,CD,Cm,rho,g,CA,ma]...
=aerodynamics(h,Maero,Vabs,rg1,Mg,alfa,beta);

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,...
    CLX,CDX,b1X,CyX,CnX]=constants;

arc = pi/180.;

% ********************************************************
% Generation of Aerodynamics of Space Vehicle
% ********************************************************

rho = rhos*exp(-bet*h);
g   = gs*(rearth/rg1)^2;

norm = rho/2.*Vabs^2*Sref;

algr = alfa/arc;
begr = beta/arc;

% ********************************************************
% X-38 Longitudinal and Lateral Aerodynamics
% ********************************************************

CL =@(alpha) -0.02075 + 0.015*alpha + 0.0254*alpha^2 + 0.0588*alpha^3 + 0.0594*alpha^4;
CD =@(alpha)  0.03553 + 0.0008*alpha - 0.0248*alpha^2 + 0.3720*alpha^3 - 0.6208*alpha^4;
Cm =@(alpha)  0.0126 - 0.0097*alpha + 0.0336*alpha^2 - 0.1431*alpha^3;
CQ =@(alpha) -0.1053*beta + 0.0265*beta^2 + 0.1008*beta^3;
Cn =@(alpha)  0.0108*beta - 0.0012*beta^2 - 0.0004*beta^3;


% CL = interp1(alX,CLX,algr,'pchip');
% CD = interp1(alX,CDX,algr,'pchip');
% Cm = interp1(alX,CmX,algr,'pchip');
% CQ = interp1(b1X,CyX,begr,'pchip');
% Cn = interp1(b1X,CnX,begr,'pchip');
Cll = 0.;

aa(1) = -CD(algr);
aa(2) =  CQ(begr);
aa(3) = -CL(algr);

CA = Maero*aa';

A = norm*Maero*aa';

ma(1) = Cll*Lref;
ma(2) = Cm(algr) *Lref;
ma(3) = Cn(begr) *Lref;

M = norm*ma';

G(1) = 0.;
G(2) = 0.;
G(3) = m*g;

MG1=Mg*G';

% ********************************************************
% End of Function Aerodynamics 
% ********************************************************
