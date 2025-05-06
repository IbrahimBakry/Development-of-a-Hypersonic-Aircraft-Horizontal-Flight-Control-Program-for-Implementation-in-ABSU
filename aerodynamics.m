% ********************************************************
% Function Aerodynamics
% Aerodynamic Forces and Moments in Body Fixed Coordinates
% ********************************************************

function [A,M,MG1,CL,CD,Cm]=aerodynamics(h,Maero,Vabs,rg1,Mg,alfa,beta)
[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,CLX,CDX,b1X,CyX,CnX]=constants;

arc = pi/180.;

% ********************************************************
% Generation of Aerodynamics of Space Vehicle
% ********************************************************

% Density Model (1)
% H=10351.8-(0.0368512).*h-(1.02368e-5).*h.^2+(2.63363e-10).*h.^3;
% rho = rhos.*exp(-h./H);
% % Density Model (2) - Old
rho = rhos*exp(-bet*h);
% %End Density Model

% Gravity Model with Hight
g = gs.*(rearth./rg1).^2;  % rg1 = h + rearth
% End Gravity Model

% Force Normalization
norm = rho/2.*Vabs^2*Sref;
% End Force Normalization

% Convert Degree-to-Radian
algr = alfa/arc;
begr = beta/arc;
% End Convert Degree-to-Radian

% ********************************************************
% Bank Angle Factor
% ********************************************************

mue = 0.4;

% ********************************************************
% Longitudinal and Lateral Aerodynamics
% ********************************************************
CL = interp1(alX,CLX,algr,'cubic')*mue;
CD = interp1(alX,CDX,algr,'cubic');
Cm = interp1(alX,CmX,algr,'cubic');
CQ = interp1(b1X,CyX,begr,'cubic');
Cn = interp1(b1X,CnX,begr,'cubic');
Cll = 0.;

% Forces Value
aa(1) = -CD;
aa(2) =  CQ;
aa(3) = -CL;

% Non-dimentional Forces
A = norm*Maero*aa';

% Moments Values
ma(1) = Cll*Lref;
ma(2) = Cm *Lref;
ma(3) = Cn *Lref;

% Non-dimentional Moments
M = norm*ma';

% Gravity Force in Inertial System coordinate
G(1) = 0.;
G(2) = 0.;
G(3) = m*g;

% Gravity Forces in Body coordinate
MG1=Mg*G';

% ********************************************************
% End of Function Aerodynamics 
% ********************************************************