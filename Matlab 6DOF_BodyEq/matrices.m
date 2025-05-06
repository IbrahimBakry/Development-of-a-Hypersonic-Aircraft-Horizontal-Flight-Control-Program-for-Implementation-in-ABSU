% ********************************************************
% Function "Matrices" Creates the Transformation Matrices%
% ********************************************************

function [Maero,Vabs,Mg,alfa,beta]...
= matrices(u1,v1,w1,phis1,thetas1,psis1);

Vabs=sqrt(u1.*u1+v1.*v1+w1.*w1);
uvb  = sqrt(u1.*u1+v1.*v1+w1.*w1);
alfa = atan(w1/uvb);
beta = atan(v1/u1);

%*********** alfa & beta limits ***********
arc = pi/180; 
if alfa <= -3*arc, alfa = -3*arc; end
if alfa >= 21*arc, alfa = 21*arc; end
if beta <= -10*arc, beta = -10*arc; end
if beta >=  10*arc, beta =  10*arc; end
%*********** End alfa & beta limits ***********


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
% End of Function Matrices
% ********************************************************