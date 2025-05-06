function [A,M,MG1 ,CL,CD,Cm] = aerodynamics(h,Maero,Vabs,rg1,Mg,alfa,beta)
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =
% Function Aerodynamics
% Provision of Aerodynamic Forces and Moments in Body Fixed Coordinates
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,a1X,CmX,CLX,CDX,b1 X,CyX,CnX] = constants;
arc = pi/180.;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Generation of Aerodynamics of Space Vehicle  (check this part, is the eq. correcect or not)
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
rho = rhos * exp(-bet*h);
g = gs * (rearth/rg1) * 2;
norm = rho/2.*Vabs * 2*Sref;
algr = alfa/arc;
begr = beta/arc;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Bank Angle Factor
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
mue = 0.4;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% X-38 Longitudinal and Lateral Aero dynamics
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
CL = interpl(a1X,CLX,algr,'cubic')* mue;
CD = interpl(a1X,CDX,algr,'cubic');
Cm = interpl(a1X,CmX,algr,'cubic');
CQ = interpl(b1X,CyX,begr,'cubic');
Cn = interpl(b1X,CnX,begr,'cubic');
Cll = 0.;
aa(1) = -CD;
aa(2) = CQ ;
aa(3) = -CL ;
A = norm * Maero * aa';
ma(1) = Cll* Lref;
ma(2) = Cm * Lref;
ma(3) = Cn * Lref;
M = norm * ma';
G(1) = 0.;
G(2) = 0.;
G(3) = m*g ;
MG1 = Mg*G';

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **
% End of Function Aerodynamics
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **