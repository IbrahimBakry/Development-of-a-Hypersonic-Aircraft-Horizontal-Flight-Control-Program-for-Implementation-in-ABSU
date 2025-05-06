% ********************************************************
% Function "Initiala" Provides and Generates Intitial Data
% ********************************************************

function [Vtot,gamma,chi,phisn,thetasn,psisn,hn,phin,thetan,mf]=initiala

% ********************************************************
% Matrix Transforming from Geodetic to Body Fixed System
% ********************************************************

arc     = pi/180.;

Vtot    = 1890;
chi    = 0*arc;
gamma   = 0*arc;
phisn   = 0.*arc;
thetasn = 10.0*arc;
psisn   = 90.*arc;
hn      = 30000;
phin    = 55*arc;
thetan  = 33*arc;
mf      = 81646.63;

% Mg(1,1) = cos(psisn)*cos(thetasn);
% Mg(1,2) = sin(psisn)*cos(thetasn);
% Mg(1,3) =-sin(thetasn);
% Mg(2,1) = sin(phisn)*sin(thetasn)*cos(psisn) - cos(phisn)*sin(psisn);
% Mg(2,2) = sin(phisn)*sin(thetasn)*sin(psisn) + cos(phisn)*cos(psisn);
% Mg(2,3) = sin(phisn)*cos(thetasn);
% Mg(3,1) = cos(phisn)*sin(thetasn)*cos(psisn) + sin(phisn)*sin(psisn);
% Mg(3,2) = cos(phisn)*sin(thetasn)*sin(psisn) - sin(phisn)*cos(psisn);
% Mg(3,3) = cos(phisn)*cos(thetasn);

% uvg     = cos(gamma)*Vtot;
% Vg(1)   = uvg*sin(chie);
% Vg(3)   = -sin(gamma)*Vtot;
% Vg(2)   = uvg*cos(chie);
% 
% Vb = Mg*Vg';
% 
% u1   = Vb(1);
% v1   = Vb(2);
% w1   = Vb(3);
% uvb = sqrt(u1.*u1 + v1.*v1 + w1.*w1);
% alfa = atan(w1/uvb);

% ********************************************************
% End of Function Initiala
% ********************************************************
