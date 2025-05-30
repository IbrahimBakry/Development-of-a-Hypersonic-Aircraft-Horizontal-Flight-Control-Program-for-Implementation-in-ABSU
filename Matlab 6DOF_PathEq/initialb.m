% ********************************************************
% Function "Initialb" Generates Intitial Data for Omega
% ********************************************************

function [omb,ome] = initialb(phis1,thetas1,psis1,h1,phi1,theta1,omegae)

[Vabs,gamma,chi,phisn,thetasn,psisn,hn,phin,thetan,mf]=initiala;
[m,Sref,Lref,bref,Ac,a,rhos,gs,bet,omegae,rearth,V00]=constants;

% ********************************************************
% Matrix Transforming from Geodetic to Body Fixed System
% ********************************************************


R = h1 + rearth;

Mg(1,1) = cos(psis1)*cos(thetas1);
Mg(1,2) = sin(psis1)*cos(thetas1);
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

h1p     =-Vabs*sin(gamma);
phi1p   = Vabs*cos(gamma)*cos(chi)/R;
theta1p = Vabs*cos(gamma)*sin(chi)/(R*cos(phi1));

om2(1)  = theta1p*cos(phi1);
om2(2)  =-phi1p;
om2(3)  =-theta1p*sin(phi1);

ome     = Mg*om2';

% ********************************************************
% End of Function Initialb
% ********************************************************
