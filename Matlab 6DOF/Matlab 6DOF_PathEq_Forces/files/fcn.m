function dyt = fcn(y,in)
%*********************************************************************
% Constants
%*********************************************************************
m0      = 136077.7;% Takeoff weight;
omegae = 0.00007292;
rearth = 6.38E6;

%*********************************************************************
% Assigning Variables
%*********************************************************************
V1      = y(1);    gamma1  = y(2);    chi1    = y(3);
h1      = y(4);   phi1    = y(5);   theta1  = y(6);
mf      = y(7);   

psiV    = y(8);    psigamma= y(9);    psichi  = y(10);
psih      = y(11); psiphi= y(12);     psitheta  = y(13);
psim      = y(14); 

Thrust = in(1);   Thrustd = in(2);    I = in(3);
F = in(4:6);      F_v = in (7:9);     g = in(10);
rho = in(10);


% % ********************************************************
% % Estimation of Position
% % ********************************************************
R = h1 + rearth;
h1p     = V1*sin(gamma1);
phi1p   = V1*cos(gamma1)*cos(chi1)/R;
theta1p = V1*cos(gamma1)*sin(chi1)/(R*cos(phi1));

% ***************************************************************
% Mass Consuming Equation 
% ***************************************************************

mfd = - Thrust/(g*I);
m = m0 - mf; 

% ***************************************************************
% Set up of Right Hand Side Eauations
% ***************************************************************

D = F(1); Y = F(2); L = F(3); mu = 0; 
BB = zeros(1,3);
BB(1)=((Thrust-D)/m)-g*sin(gamma1)+omegae^2*R*cos(phi1)^2*...
    (sin(gamma1)-cos(gamma1)*tan(phi1)*sin(chi1));
BB(2)=(1/V1)*((L*cos(mu)-Y*sin(mu))/m - g*cos(gamma1)+(V1^2*...
    cos(gamma1))/R+2*omegae*V1*cos(phi1)*cos(chi1)+omegae^2*R*...
    cos(phi1)^2*(cos(gamma1)+sin(gamma1)*tan(phi1)*sin(chi1)));
BB(3)=(1/V1)*((L*sin(mu)+Y*cos(mu))/(m*cos(gamma1))-(V1^2/R)*cos(gamma1)...
    *cos(chi1)*tan(phi1)+2*omegae*V1*(tan(gamma1)*cos(phi1)*sin(chi1)...
    -sin(phi1))-(omegae^2*R/cos(gamma1))*sin(phi1)*cos(phi1)*cos(chi1));

x1 = [BB(1) BB(2) BB(3)]';

% ***************************************************************
% Complete Set of the 7 Ordinary Differential Equations
% ***************************************************************
dydt = zeros(1,14);
dydt(1) = x1(1);
dydt(2) = x1(2);
dydt(3) = x1(3);

dydt(4) = h1p;
dydt(5) = phi1p;
dydt(6) = theta1p;  
dydt(7) = mfd;

%****************************************************************
% Alpha Costate 7 Differential Equations  
%**************************************************************** 
Dv = F_v(1); Yv = F_v(2); Lv = F_v(3);

psiVd = -psiV*(1/m)*(Thrustd-Dv) - psigamma*((1/(m*V1))*Lv +...
    cos(gamma1)/R + (1/V1^2)*(g*cos(gamma1)-R*omegae^2*cos(phi1)^2*...
    (sin(gamma1)*tan(phi1)*sin(chi1)+cos(gamma1))-L/m))-psichi*...
    (-cos(gamma1)*tan(phi1)*cos(chi1)/R + (omegae^2*R*sin(phi1)*...
    cos(phi1)*cos(chi1))/(cos(gamma1)*V1^2))...
    - psih*sin(gamma1) - psiphi*...
    (cos(gamma1)*cos(chi1)/R) - psitheta*(cos(gamma1)*sin(chi1)/...
    (R*cos(phi1)));

psigammad = -psiV*(-g*cos(gamma1)+R*omegae^2*cos(phi1)^2*(sin(gamma1)*...
    tan(phi1)*sin(chi1)+cos(gamma1))) - psigamma*(R*omegae^2*cos(phi1)...
    ^2*(cos(gamma1)*tan(phi1)*sin(chi1)-sin(gamma1))/V1+sin(gamma1)*...
    (g/V1-V1/R)) - psichi*(-R*omegae^2*sin(phi1)*cos(phi1)*cos(chi1)*...
    sin(gamma1)/(V1*cos(gamma1)^2)+2*omegae*cos(phi1)*sin(chi1)*...
    (tan(gamma1)^2+1)+V1*sin(gamma1)*cos(chi1)*tan(phi1)/R) - psih*...
    (V1*cos(gamma1)) + psiphi*...
    (V1*sin(gamma1)*cos(chi1)/R) + psitheta*(V1*sin(gamma1)*sin(chi1)/...
    (R*cos(phi1)));

psichid = psiV*(R*omegae^2*cos(phi1)^2*cos(gamma1)*tan(phi1)*cos(chi1))-...
    psigamma*(-2*omegae*cos(phi1)*sin(chi1)+R*omegae^2*cos(phi1)^2*...
    sin(gamma1)*tan(phi1)*cos(chi1)/V1) - psichi*(V1*cos(gamma1)*...
    sin(chi1)*tan(phi1)/R+2*omegae*tan(gamma1)*cos(phi1)*cos(chi1)+...
    omegae^2*R*sin(phi1)*cos(phi1)*sin(chi1)/(V1*cos(gamma1))) +...
    psiphi*(V1*cos(gamma1)*...
    sin(chi1)/R) - psitheta*(V1*cos(gamma1)*cos(chi1)/(R*cos(phi1)));

psihd = -psiV*(omegae^2*cos(phi1)^2*(-cos(gamma1)*tan(phi1)*sin(chi1)+...
    sin(gamma1))) - psigamma*(-V1*cos(gamma1)/R^2+omegae^2*cos(phi1)^2*(...
    cos(gamma1)+sin(gamma1)*tan(phi1)*sin(chi1))) - psichi*(V1*cos(...
    gamma1)*cos(chi1)*tan(phi1)/R^2-omegae^2*sin(phi1)*cos(phi1)*...
    cos(chi1)/(V1*cos(gamma1))) + psiphi*(V1*cos(gamma1)*cos(chi1)/R^2)...
    + psitheta*(V1*cos(gamma1)*cos(chi1)/(R^2*cos(phi1)));

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

dydt(8) = psiVd;
dydt(9) = psigammad;
dydt(10) = psichid;
dydt(11) = psihd;
dydt(12) = psiphid;
dydt(13) = psithetad;
dydt(14) = psimd;

dyt = double(dydt');

% ***************************************************************
% End of Sub Program RHS
% ***************************************************************