function [A,M,MG1,CA,ma,pb,qb,rb,centripetal,pp1,qq1,rr1,omb]=repeatRHS...
    (h1,Vabs,rg1,alfa,beta,omegae,phi1,u1,v1,w1,thetas1,psis1,phis1,p1,q1,r1)

% ********************************************************
% Matrix Transforming from Aerodynamic to Body Fixed System
% ********************************************************

Maero(1,1) = cos(alfa)*cos(beta);
Maero(1,2) =-cos(alfa)*sin(beta);
Maero(1,3) =-sin(alfa);
Maero(2,1) = sin(beta);
Maero(2,2) = cos(beta);
Maero(2,3) = 0.;
Maero(3,1) = sin(alfa)*cos(beta);
Maero(3,2) =-sin(alfa)*sin(beta);
Maero(3,3) = cos(alfa);

% ********************************************************
% Matrix Transforming from Geodetic to Body Fixed System
% ********************************************************

Mg(1,1) = cos(psis1)*cos(thetas1);
Mg(1,2) = sin(psis1)*cos(thetas1);
Mg(1,3) =-sin(thetas1);
Mg(2,1) = sin(phis1)*sin(thetas1)*cos(psis1) - cos(phis1)*sin(psis1);
Mg(2,2) = sin(phis1)*sin(thetas1)*sin(psis1) + cos(phis1)*cos(psis1);
Mg(2,3) = sin(phis1)*cos(thetas1);
Mg(3,1) = cos(phis1)*sin(thetas1)*cos(psis1) + sin(phis1)*sin(psis1);
Mg(3,2) = cos(phis1)*sin(thetas1)*sin(psis1) - sin(phis1)*cos(psis1);
Mg(3,3) = cos(phis1)*cos(thetas1);

% ********************************************************
% Establishment of Elimination Matrix AA
% ********************************************************

om1(1)  =  cos(phi1)*omegae;
om1(2)  =  0.;
om1(3)  = -sin(phi1)*omegae;

omb     = Mg*om1';

pb  = omb(1); qb  = omb(2); rb  = omb(3);

ommat= zeros(3,3);
ommat(1,2) = -rb; ommat(1,3) =  qb; ommat(2,1) =  rb;
ommat(2,3) = -pb; ommat(3,1) = -qb; ommat(3,2) =  pb;
ommat2 = ommat*ommat;

rg1v(1) = 0.; rg1v(2) = 0.; rg1v(3) = rg1; rg1vb = Mg*rg1v';
centripetal = -ommat2*rg1vb;

V1(1)  = u1; V1(2)  = v1; V1(3)  = w1;
V1g    = Mg'*V1';

gamma  =-asin(V1g(3)/Vabs);
chi    = atan(V1g(2)/V1g(1));

if V1g(1) < 0.
    chi = pi + atan(V1g(2)/V1g(1));
end

%************* Gamma & Chi Limits **************
arc = pi/180;
if gamma <= -45*arc;  gamma = -45*arc; end
if gamma <=  45*arc;  gamma =  45*arc; end
if chi <= -45*arc;  chi = -45*arc; end
if chi <=  45*arc;  chi =  45*arc; end
%************ End Gamma & Chi Limits **************

h1p     = Vabs*sin(gamma);
phi1p   = Vabs*cos(gamma)*cos(chi)/rg1;
theta1p = Vabs*cos(gamma)*sin(chi)/(rg1*cos(phi1));

om2(1)  = theta1p*cos(phi1);
om2(2)  =-phi1p;
om2(3)  =-theta1p*sin(phi1);

ome     = Mg*om2';
pe  = ome(1); qe  = ome(2); re  = ome(3);

pp1  = p1 - pb - pe;
qq1  = q1 - qb - qe;
rr1  = r1 - rb - re;

[A,M,MG1,CL,CD,Cm,rho,g,CA,ma]...
    =aerodynamics(h1,Maero,Vabs,rg1,Mg,alfa,beta);
