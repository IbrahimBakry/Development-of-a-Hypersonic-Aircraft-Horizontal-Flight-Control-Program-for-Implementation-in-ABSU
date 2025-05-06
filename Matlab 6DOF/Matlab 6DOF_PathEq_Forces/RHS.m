% ********************************************************
% Sub Program Providing the Right Hand Sides of the Differential 
% Equation dy/dt= by Gauss Elimination of A y(t) = B.
% ********************************************************

function dyt=RHS(t,y)

[m0,Sref,Lref,bref,Ac,a,rhos,gs,bet,omegae,rearth,V00]=constants;

[T] = inertia;

V1      = y(1);    gamma1  = y(2);    chi1    = y(3);
p1      = y(4);    q1      = y(5);    r1      = y(6);
phis1   = y(7);    thetas1 = y(8);    psis1   = y(9);
h1      = y(10);   phi1    = y(11);   theta1  = y(12);
mf      = y(13);   psiV    = y(14);   psigamma= y(15);
psichi  = y(16);   psiqt   = y(17);   psithetae = y(18);
psih      = y(19); psiphi    = y(20); psitheta  = y(21);
psim      = y(22); 



%******************************************************
% Controller
%******************************************************

kp  = 0.1;
hd  = 30000; 
C   = -kp*(h1-hd);
% % V11 = C;
h1p = C;



%******************************************************
%  Limits of Angles and Hight
%******************************************************

arc = pi/180;
if V1 <= 1525;   V1 = 1525;  end
if V1 <=  2135;  V1 =  2135; end
if gamma1 <= -45*arc;  gamma1 = -45*arc; end
if gamma1 <=  45*arc;  gamma1 =  45*arc; end
if chi1 <= -45*arc;  chi1 = -45*arc; end
if chi1 <=  45*arc;  chi1 =  45*arc; end
if h1 <= 15000, h1 = 15000; end
if h1 >= 40000, h1 = 40000; end
if phi1 <= -89*arc, phi1 = -89*arc; end
if phi1 >=  89*arc, phi1 =  89*arc; end
if theta1 <= -180*arc, theta1 = -180*arc; end
if theta1 >=  180*arc, theta1 =  180*arc; end
%******* End Limits of Angles and Hight **********

%******************************************************
% condition of alfa maximum
%******************************************************

if psiqt <= 0, psiqt = - psiqt; end
if psigamma >= 0, psigamma = - psigamma; end
if psiV <= 0, psiV = - psiV; end



[alfa] = alpha_equation(V1,psiV,psigamma,psiqt);
% [beta] = beta_equation(V1,gamma1,psichi2,psipt2,psirt2);
beta=0;


[Maero,Mg] = matrices(phis1,thetas1,psis1,alfa,beta);

R = h1 + rearth;

[F,F_v,M,M_v,MG1,rho,g,ma] = aerodynamics(h1,V1,R,Mg,alfa,beta);

% % ********************************************************
% % Establishment of Elimination Matrix AA
% % ********************************************************

% h1p     = V1*sin(gamma1);
phi1p   = V1*cos(gamma1)*cos(chi1)/R;
theta1p = V1*cos(gamma1)*sin(chi1)/(R*cos(phi1));


%***************************************************************
% Thrust Equations
%***************************************************************

Mach = V1/a;
[I,dI1,dI2,p,pd]=Engine_CT(Mach);

CTP  = polyval(p,Mach);
CTPd = polyval(pd,Mach);

Thrust  = 0.029*rho*g*V1*I*Ac*CTP;
Thrustd = 0.029*rho*g*dI1*Ac*CTP*V1 + 0.029*rho*g*I*Ac*CTPd*V1 +... 
          0.029*rho*g*(1/a)*I*Ac*CTP;

% ***************************************************************
% Mass Consuming Equation 
% ***************************************************************

mfd = - Thrust/(g*I);
m = m0 - mf; 


% ***************************************************************
% Set up of Right Hand Side (The Inhomogeneous part of the ODE)
% ***************************************************************

D = F(1); Y = F(2); L = F(3); mu = 0; 

BB(1)=((Thrust-D)/m)-g*sin(gamma1)+omegae^2*R*cos(phi1)^2*...
    (sin(gamma1)-cos(gamma1)*tan(phi1)*sin(chi1));
BB(2)=(1/V1)*((L*cos(mu)-Y*sin(mu))/m - g*cos(gamma1)+(V1^2*...
    cos(gamma1))/R+2*omegae*V1*cos(phi1)*cos(chi1)+omegae^2*R*...
    cos(phi1)^2*(cos(gamma1)+sin(gamma1)*tan(phi1)*sin(chi1)));
BB(3)=(1/V1)*((L*sin(mu)+Y*cos(mu))/(m*cos(gamma1))-(V1^2/R)*cos(gamma1)...
    *cos(chi1)*tan(phi1)+2*omegae*V1*(tan(gamma1)*cos(phi1)*sin(chi1)...
    -sin(phi1))-(omegae^2*R/cos(gamma1))*sin(phi1)*cos(phi1)*cos(chi1));

x1 = [BB(1) BB(2) BB(3)]';
x2 = [0 0 0]';

% ***************************************************************
% Complete Set of the 13 Ordinary Differential Equations
% ***************************************************************

dydt(1) = x1(1);
dydt(2) = x1(2);
dydt(3) = x1(3);
dydt(4) = 0;
dydt(5) = 0;
dydt(6) = 0;

dydt(7)  = 0;
dydt(8)  = 0;
dydt(9)  = 0;
dydt(10) = h1p;
dydt(11) = phi1p;
dydt(12) = theta1p;  
dydt(13) = mfd;

%****************************************************************
% Alpha Costate 9 Differential Equations  
%**************************************************************** 
Dv = F_v(1); Yv = F_v(2); Lv = F_v(3);
Mv = M_v(2);  Iy = T(2,2);

psiVd = -psiV*(1/m)*(Thrustd-Dv) - psigamma*((1/(m*V1))*Lv +...
    cos(gamma1)/R + (1/V1^2)*(g*cos(gamma1)-R*omegae^2*cos(phi1)^2*...
    (sin(gamma1)*tan(phi1)*sin(chi1)+cos(gamma1))-L/m))-psichi*...
    (-cos(gamma1)*tan(phi1)*cos(chi1)/R + (omegae^2*R*sin(phi1)*...
    cos(phi1)*cos(chi1))/(cos(gamma1)*V1^2)) - psiqt*(Mv/Iy) -...
    psithetae*(cos(gamma1)*cos(chi1)/R) - psih*sin(gamma1) - psiphi*...
    (cos(gamma1)*cos(chi1)/R) - psitheta*(cos(gamma1)*sin(chi1)/...
    (R*cos(phi1)));

psigammad = -psiV*(-g*cos(gamma1)+R*omegae^2*cos(phi1)^2*(sin(gamma1)*...
    tan(phi1)*sin(chi1)+cos(gamma1))) - psigamma*(R*omegae^2*cos(phi1)...
    ^2*(cos(gamma1)*tan(phi1)*sin(chi1)-sin(gamma1))/V1+sin(gamma1)*...
    (g/V1-V1/R)) - psichi*(-R*omegae^2*sin(phi1)*cos(phi1)*cos(chi1)*...
    sin(gamma1)/(V1*cos(gamma1)^2)+2*omegae*cos(phi1)*sin(chi1)*...
    (tan(gamma1)^2+1)+V1*sin(gamma1)*cos(chi1)*tan(phi1)/R) + psithetae...
    *(V1*sin(gamma1)*cos(chi1)/R) - psih*(V1*cos(gamma1)) + psiphi*...
    (V1*sin(gamma1)*cos(chi1)/R) + psitheta*(V1*sin(gamma1)*sin(chi1)/...
    (R*cos(phi1)));

psichid = psiV*(R*omegae^2*cos(phi1)^2*cos(gamma1)*tan(phi1)*cos(chi1))-...
    psigamma*(-2*omegae*cos(phi1)*sin(chi1)+R*omegae^2*cos(phi1)^2*...
    sin(gamma1)*tan(phi1)*cos(chi1)/V1) - psichi*(V1*cos(gamma1)*...
    sin(chi1)*tan(phi1)/R+2*omegae*tan(gamma1)*cos(phi1)*cos(chi1)+...
    omegae^2*R*sin(phi1)*cos(phi1)*sin(chi1)/(V1*cos(gamma1))) +...
    psithetae*(V1*cos(gamma1)*sin(chi1)/R) + psiphi*(V1*cos(gamma1)*...
    sin(chi1)/R) - psitheta*(V1*cos(gamma1)*cos(chi1)/(R*cos(phi1)));

psiqrd = 0;

psithetaed = 0;

psihd = -psiV*(omegae^2*cos(phi1)^2*(-cos(gamma1)*tan(phi1)*sin(chi1)+...
    sin(gamma1))) - psigamma*(-V1*cos(gamma1)/R^2+omegae^2*cos(phi1)^2*(...
    cos(gamma1)+sin(gamma1)*tan(phi1)*sin(chi1))) - psichi*(V1*cos(...
    gamma1)*cos(chi1)*tan(phi1)/R^2-omegae^2*sin(phi1)*cos(phi1)*...
    cos(chi1)/(V1*cos(gamma1))) + psithetae*(V1*cos(gamma1)*cos(chi1)...
    /R^2) + psiphi*(V1*cos(gamma1)*cos(chi1)/R^2) + psitheta*(V1*...
    cos(gamma1)*cos(chi1)/(R^2*cos(phi1)));

psiphid = -psiV*R*omegae^2*(-cos(gamma1)*sin(phi1)^2*sin(chi1)+2*...
    cos(gamma1)*sin(phi1)^2*sin(chi1)-cos(gamma1)*cos(phi1)^2*sin(chi1)-...
    sin(gamma1)*sin(2*phi1)) - psigamma*(-2*omegae*sin(phi1)*cos(chi1)...
    +(R*omegae^2/V1)*(-sin(2*phi1)*cos(gamma1)-sin(2*phi1)*sin(gamma1)*...
    tan(phi1)*sin(chi1)+cos(phi1)^2*sin(gamma1)*sin(chi1)+sin(phi1)^2*...
    sin(gamma1)*sin(chi1))) - psichi*(-V1*cos(gamma1)*cos(chi1)*(tan...
    (phi1)^2+1)/R-2*omegae*(tan(gamma1)*sin(phi1)*sin(chi1)+cos(phi1))...
    +R*omegae^2*cos(chi1)*(sin(phi1)^2-cos(phi1)^2)/(V1*cos(gamma1))) +...
    psitheta*(V1*cos(gamma1)*sin(chi1)*sin(phi1)/(R*cos(phi1)^2));

psithetad = 0;

psimd = - psiV*((-Thrust+D)/m^2) + psigamma*(L/(V1*m^2));

dydt(14) = psiVd;
dydt(15) = psigammad;
dydt(16) = psichid;
dydt(17) = psiqrd;
dydt(18) = psithetaed;
dydt(19) = psihd;
dydt(20) = psiphid;
dydt(21) = psithetad;
dydt(22) = psimd;

% %****************************************************************
% % Beta Costate 10 Differential Equations  
% %****************************************************************
% Ix = T(1,1); Iz = T(3,3); Ixz = T(1,3); II = (Ix*Iz-Ixz^2);
% Lmv = M_v(1); Nv = M_v(3); 
% psiV2d = - psiV2*(1/m)*(Thrustd-Dv) - psigamma2*(+...
%     cos(gamma1)/R + (1/V1^2)*(g*cos(gamma1)-R*omegae^2*cos(phi1)^2*...
%     (sin(gamma1)*tan(phi1)*sin(chi1)+cos(gamma1)))) - psichi2*...
%     (Yv/(m*cos(gamma1)*V1)-cos(gamma1)*tan(phi1)*cos(chi1)/R + (1/V1^2)...
%     *(-Y/(m*cos(gamma1))+omegae^2*R*sin(phi1)*cos(phi1)*cos(chi1)/...
%     (cos(gamma1)))) - psipt2*(Lmv*Iz/II + Nv*Ixz/II) - psirt2*(Nv*II...
%     +Ixz*(Lmv*Iz+Nv*Ixz))/(Iz*II) - psipsie2*(cos(gamma1)*sin(chi1)...
%     *sin(phi1)/(R*cos(phi1))) + psiphie2*(cos(psis1)*cos(gamma1)*sin(chi1)...
%     /R) - psiphi2*(cos(gamma1)*cos(chi1)/R) - psitheta2*(cos(gamma1)*...
%     sin(chi1)/(R*cos(phi1)));
% 
% psigamma2d = -psiV2*(-g*cos(gamma1)+R*omegae^2*cos(phi1)^2*(sin(gamma1)*...
%     tan(phi1)*sin(chi1)+cos(gamma1))) - psigamma2*(R*omegae^2*cos(phi1)...
%     ^2*(cos(gamma1)*tan(phi1)*sin(chi1)-sin(gamma1))/V1+sin(gamma1)*...
%     (g/V1-V1/R)) - psichi2*(-R*omegae^2*sin(phi1)*cos(phi1)*cos(chi1)*...
%     sin(gamma1)/(V1*cos(gamma1)^2)+2*omegae*cos(phi1)*sin(chi1)*...
%     (tan(gamma1)^2+1)+V1*sin(gamma1)*cos(chi1)*tan(phi1)/R + sin(gamma1)...
%     *Y/(V1*m*cos(gamma1)^2)) + psipsie2*(V1*sin(gamma1)*cos(chi1)*...
%     sin(phi1)/(R*cos(phi1))) - psiphie2*(V1*cos(phis1)*sin(gamma1)*...
%     sin(chi1)/R) + psiphi2*(V1*sin(gamma1)*cos(chi1)/R) + psitheta2*...
%     (V1*sin(gamma1)*sin(chi1)/(R*cos(phi1)));
% 
% psichi2d = psiV2*(R*omegae^2*cos(phi1)^2*cos(gamma1)*tan(phi1)*cos(chi1...
%     )) - psigamma2*(-2*omegae*cos(phi1)*sin(chi1)+R*omegae^2*cos(phi1)^2*...
%     sin(gamma1)*tan(phi1)*cos(chi1)/V1) - psichi2*(V1*cos(gamma1)*...
%     sin(chi1)*tan(phi1)/R+2*omegae*tan(gamma1)*cos(phi1)*cos(chi1)+...
%     omegae^2*R*sin(phi1)*cos(phi1)*sin(chi1)/(V1*cos(gamma1))) -...
%     psipsie2*(V1*cos(gamma1)*cos(chi1)*sin(phi1)/(R*cos(phi1)))...
%     + psiphie2*(V1*cos(psis1)*cos(gamma1)*cos(chi1)/R) + psiphi2*(V1*...
%     cos(gamma1)*sin(chi1)/R) - psitheta2*(V1*cos(gamma1)*cos(chi1)/(R*...
%     cos(phi1)));
% 
% psipt2d = - psiphie2;
% 
% psirt2d = - psipsie2*cos(phis1);
% 
% psiphie2d = - psipsie2*(q1*cos(phis1)-r1*sin(phis1));
% 
% psipsie2d = - psiphie2*(omegae*sin(psis1)*cos(phi1)+V1*sin(psis1)*...
%     cos(gamma1)*sin(chi1)/R);
% 
% psiphi2d = -psiV2*R*omegae^2*(-cos(gamma1)*sin(phi1)^2*sin(chi1)+2*...
%     cos(gamma1)*sin(phi1)^2*sin(chi1)-cos(gamma1)*cos(phi1)^2*sin(chi1)-...
%     sin(gamma1)*sin(2*phi1)) - psigamma2*(-2*omegae*sin(phi1)*cos(chi1)...
%     +(R*omegae^2/V1)*(-sin(2*phi1)*cos(gamma1)-sin(2*phi1)*sin(gamma1)*...
%     tan(phi1)*sin(chi1)+cos(phi1)^2*sin(gamma1)*sin(chi1)+sin(phi1)^2*...
%     sin(gamma1)*sin(chi1))) - psichi2*(-V1*cos(gamma1)*cos(chi1)*(tan...
%     (phi1)^2+1)/R-2*omegae*(tan(gamma1)*sin(phi1)*sin(chi1)+cos(phi1))...
%     +R*omegae^2*cos(chi1)*(sin(phi1)^2-cos(phi1)^2)/(V1*cos(gamma1)))-...
%     psipsie2*(omegae*cos(phi1)+V1*cos(gamma1)*sin(chi1)*sin(phi1)^2/...
%     (R*cos(phi1)^2)+V1*cos(gamma1)*sin(chi1)/R) - psiphie2*(omegae*...
%     cos(psis1)*sin(phi1)) - psitheta2*(V1*cos(gamma1)*sin(chi1)*...
%     sin(phi1)/(R*cos(phi1)^2));
% 
% psitheta2d = 0;
% 
% psim2d = -psiV2*((1/m^2)*(-Thrust+D)) + psichi2*(Y/(V1*m^2*cos(gamma1)));
% 
% dydt(23) = psiV2d;
% dydt(24) = psigamma2d;
% dydt(25) = psichi2d;
% dydt(26) = psipt2d;
% dydt(27) = psirt2d;
% dydt(28) = psiphie2d;
% dydt(29) = psipsie2d;
% dydt(30) = psiphi2d;
% dydt(31) = psitheta2d;
% dydt(32) = psim2d;
dyt = double(dydt');

t
% ***************************************************************
% End of Sub Program RHS
% ***************************************************************