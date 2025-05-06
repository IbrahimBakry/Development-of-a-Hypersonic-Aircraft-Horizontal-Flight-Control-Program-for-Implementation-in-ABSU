% ********************************************************
% Function "Initialb" Generates Intitial Data for Omega
% ********************************************************
function [omb,ome]=initialb(u1,v1,w1,phis1,thetas1,psis1,h1,phi1,theta1,omegae)
[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,CLX,CDX,b1X,CyX,CnX]=constants;

% ********************************************************
% Matrix Transforming from Geodetic to Body Fixed System
% ********************************************************

Vabs=sqrt(u1.*u1+v1.*v1+w1.*w1);
rg1 = h1 + rearth;

Mg(1,1) = cos(psis1)*cos(thetas1);   % phis1 = phi', thetas1 = theta',
Mg(1,2) = sin(psis1)*cos(thetas1);   % psis1 = psi'
Mg(1,3) =-sin(thetas1);
Mg(2,1) = sin(phis1)*sin(thetas1)*cos(psis1) - cos(phis1)*sin(psis1);
Mg(2,2) = sin(phis1)*sin(thetas1)*sin(psis1) + cos(phis1)*cos(psis1);
Mg(2,3) = sin(phis1)*cos(thetas1);
Mg(3,1) = cos(phis1)*sin(thetas1)*cos(psis1) + sin(phis1)*sin(psis1);
Mg(3,2) = cos(phis1)*sin(thetas1)*sin(psis1) - sin(phis1)*cos(psis1);
Mg(3,3) = cos(phis1)*cos(thetas1);

om1(1)  =  cos(phi1)*omegae;
om1(2)  =  0.;
om1(3)  = -sin(phi1)*omegae;

omb     = Mg*om1';

V1(1)  = u1;
V1(2)  = v1;
V1(3)  = w1;

V1g    = Mg'*V1';

gamma  =-asin(V1g(3)/Vabs);
chi    = atan(V1g(2)/V1g(1));

if V1g(1) < 0.
    chi = pi + atan(V1g(2)/V1g(1));
end

h1p     =-Vabs*sin(gamma);
phi1p   = Vabs*cos(gamma)*cos(chi)/rg1;
theta1p = Vabs*cos(gamma)*sin(chi)/(rg1*cos(phi1));

om2(1)  = theta1p*cos(phi1);
om2(2)  =-phi1p;
om2(3)  =-theta1p*sin(phi1); % phi1 = phi

ome     = Mg*om2';

% ********************************************************
% End of Function Initialb
% ********************************************************