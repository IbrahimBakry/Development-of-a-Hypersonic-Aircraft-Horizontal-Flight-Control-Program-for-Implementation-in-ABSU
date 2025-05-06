% ********************************************************
% Sub Program Providing the Right Hand Sides of the Differential 
% Equation dy/dt= by Gauss Elimination of A y(t) = B.
% ********************************************************

function dydt=RHS(t,y)

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,...
    CLX,CDX,b1X,CyX,CnX]=constants;

[T] = inertia;

u1      = y(1);
v1      = y(2);
w1      = y(3);
p1      = y(4);
q1      = y(5);
r1      = y(6);
phis1   = y(7);
thetas1 = y(8);
psis1   = y(9);
h1      = y(10);
phi1    = y(11);
theta1  = y(12);
psiu =      y(13);
psiw =      y(14);
psiqt =     y(15);
psithetae = y(16);
psih =      y(17);
psiphi =    y(18);
psim =      y(19);

%******** Limits of Angles and Hight ************
arc = pi/180;
if phis1 <= -45*arc, phis1 = -45*arc; end
if phis1 >=  45*arc, phis1 =  45*arc; end
if thetas1 <= -45*arc, thetas1 = -45*arc; end
if thetas1 >=  45*arc, thetas1 =  45*arc; end
if psis1 <= -180*arc, psis1 = -180*arc; end
if psis1 >=  180*arc, psis1 =  180*arc; end
if h1 <= 20000, h1 = 20000; end
if h1 >= 40000, h1 = 40000; end
if phi1 <= -89*arc, phi1 = -89*arc; end
if phi1 >=  89*arc, phi1 =  89*arc; end
if theta1 <= -180*arc, theta1 = -180*arc; end
if theta1 >=  180*arc, theta1 =  180*arc; end
%******* End Limits of Angles and Hight **********

[Maero,Vabs,Mg,alfa,beta]=matrices(u1,v1,w1,phis1,thetas1,psis1);

rg1 = h1 + rearth;

[A,M,MG1,CL,CD,Cm,rho,g,CA,ma]...
    =aerodynamics(h1,Maero,Vabs,rg1,Mg,alfa,beta);

% ********************************************************
% Establishment of Elimination Matrix AA
% ********************************************************

AA=zeros(6,6);
AA(1,1)=   m; AA(2,2)=   m; AA(3,3)=   m;
AA(4,4)=   T(1,1);                   AA(4,6)=   T(1,3);
                   AA(5,5)=   T(2,2);
AA(6,4)=   T(3,1);                   AA(6,6)=   T(3,3);

om1(1)  =  cos(phi1)*omegae;
om1(2)  =  0.;
om1(3)  = -sin(phi1)*omegae;

omb     = Mg*om1';

pb  = omb(1);
qb  = omb(2);
rb  = omb(3);

ommat= zeros(3,3);
                   ommat(1,2) = -rb; ommat(1,3) =  qb;
ommat(2,1) =  rb;                    ommat(2,3) = -pb;
ommat(3,1) = -qb;  ommat(3,2) =  pb;

ommat2 = ommat*ommat;

rg1v(1) = 0.;
rg1v(2) = 0.;
rg1v(3) = rg1;

rg1vb = Mg*rg1v';

centripetal = -ommat2*rg1vb;

V1(1)  = u1;
V1(2)  = v1;
V1(3)  = w1;

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

pe  = ome(1);
qe  = ome(2);
re  = ome(3);

pp1  = p1 - pb - pe;
qq1  = q1 - qb - qe;
rr1  = r1 - rb - re;

%***************************************************************
% Thrust Equations
%***************************************************************

Thrust = 10*9.81*3146.15;

%***************************************************************
% Solving Alpha Equation of Optimization
%***************************************************************

Cx = CA(1); Cy = CA(2); Cz = CA(3);
Cl = ma(1); Cm = ma(2); Cn = ma(3);


% alfa = double(alpha_equation(Vabs,rho,g,gamma,rg1,phi1,psiu,psiw,psiqt,alfa,Maero,beta));
% 
% % alfa limits
% arc = pi/180; 
% if alfa <= -3*arc, alfa = -3*arc; end
% if alfa >= 21*arc, alfa = 21*arc; end

% % repeating stage 1
% [A,M,MG1,CA,ma,pb,qb,rb,centripetal,pp1,qq1,rr1,omb]=repeatRHS...
%     (h1,Vabs,rg1,alfa(i+1),beta,omegae,phi1,u1,v1,w1,thetas1,psis1,phis1,p1,q1,r1);

% ***************************************************************
% Set up of Right Hand Side (The Inhomogeneous part of the ODE)
% ***************************************************************

BB(1)=A(1)+Thrust+MG1(1)-m*(q1*w1-r1*v1)-m*(qb*w1-rb*v1) - m*centripetal(1);
BB(2)=A(2)+MG1(2)-m*(r1*u1-p1*w1)-m*(rb*u1-pb*w1) - m*centripetal(2);
BB(3)=A(3)+MG1(3)-m*(p1*v1-q1*u1)-m*(pb*v1-qb*u1) - m*centripetal(3);
BB(4)=M(1)-T(1,3)*p1*q1+(T(2,2)-T(3,3))*q1*r1;
BB(5)=M(2)-T(1,3)*(r1^2-p1^2)+ (T(3,3)-T(1,1))*r1*p1;
BB(6)=M(3)-T(1,3)*q1*r1+(T(1,1)-T(2,2))*p1*q1;

X=inv(AA)*BB';

% ***************************************************************
% Rate of Change of the Euler Angles Describing the Attitude 
% of the Space Vehicle
% ***************************************************************

phis1p   = pp1+sin(phis1)*tan(thetas1)*qq1+cos(phis1)*tan(thetas1)*rr1;
thetas1p = cos(phis1)*qq1-sin(phis1)*rr1;
psis1p   = sin(phis1)/cos(thetas1)*qq1+cos(phis1)/cos(thetas1)*rr1;

% ***************************************************************
% Complete Set of the 12 Ordinary Differential Equations
% ***************************************************************

dydt=X;

dydt(7)  = phis1p;
dydt(8)  = thetas1p;
dydt(9)  = psis1p;
dydt(10) = h1p;
dydt(11) = phi1p;
dydt(12) = theta1p;  

%**************************************************************
% Parameters of Optimization 
%***************************************************************

Iy = T(2,2); qt = qq1;
pe = omb(1); qe = omb(2); re = omb(3);
R = rg1; u = V1(1); v = V1(2); w = V1(3); V = Vabs;
phi = phi1; thetae = thetas1; h = h1;

Tx = Thrust;

rho = rhos*exp(-bet*h1);
g   = gs*(rearth/rg1)^2;


%****************************************************************
% Alpha Costate 7 Differential Equations  
%****************************************************************

psiud = -psiu*(rho*u*Sref*Cx/m)-psiw*(rho*u*Sref*Cz/m+qt)...
    -psiqt*(rho*u*Sref*Lref*Cm/Iy)-(1/R)*(psithetae+psiphi);

psiwd = -psiu*(rho*w*Sref*Cx/m-qt)-psiw*(rho*w*Sref*Cz/m)...
    -psiqt*(rho*w*Sref*Lref*Cm/Iy)+psih;

psiqtd = psiu*w - psiw*u - psithetae;

psithetaed = psiu*(g*cos(thetae)+R*omegae^2*cos(phi)...
    *(cos(thetae)*cos(phi)+sin(thetae)*sin(phi))) - psiw*(-g*sin(thetae)...
    +R*omegae^2*cos(phi)*(cos(thetae)*sin(phi)-sin(thetae)*cos(phi)));

psihd = psiu*(re^2*sin(thetae)+re*pe*cos(thetae)) - psiw*(re*pe...
    *sin(thetae)+pe^2*cos(thetae)) + (V*cos(gamma)/R^2)*(psithetae+psiphi);
    
psiphid = -psiu*(R*omegae^2*(2*cos(thetae)*cos(phi)^2+2*cos(phi)...
    *sin(thetae)*sin(phi)-cos(thetae))) - psiw*(R*omegae^2*(2*sin(thetae...
    )*cos(phi)^2-2*cos(phi)*cos(thetae)*sin(phi)-sin(thetae)));

psimd = psiu*((0.5*rho*V^2*Sref*Cx+Tx)/m^2) + ...
    psiw*(0.5*rho*V^2*Sref*Cz/m^2);


dydt(13) = psiud;
dydt(14) = psiwd;
dydt(15) = psiqtd;
dydt(16) = psithetaed;
dydt(17) = psihd;
dydt(18) = psiphid;
dydt(19) = psimd;

t
% ***************************************************************
% End of Sub Program RHS
% ***************************************************************