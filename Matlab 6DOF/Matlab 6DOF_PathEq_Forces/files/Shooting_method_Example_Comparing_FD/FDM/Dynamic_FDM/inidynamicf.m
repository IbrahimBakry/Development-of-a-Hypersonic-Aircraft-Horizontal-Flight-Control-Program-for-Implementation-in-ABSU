function [En,E0,EM] = inidynamicf(k,dt,y)
%*********************************************************************
% Constants
%*********************************************************************
m0     = 136077.7;% Takeoff weight;
omegae = 0.00007292;
rearth = 6.38E6;
Sref = 21.672;
rhos = 1.225;
gs = 9.80665;
Ac = 27.87;
a  = 305;
mu = 0;

%*********************************************************************
% Assigning Variables
%*********************************************************************
V2 = y(k,1); V1 = 0;
gamma2 = y(k,2); gamma1 = 0;
chi2 = y(k,3); chi1 = 0;
h2 = y(k,4); h1 = 0;
phi2 = y(k,5); phi1 = 0;
theta2 = y(k,6); theta1 = 0;
mf2 = y(k,7); mf1 = 0;

%*********************************************************************
% Mid-Value
%*********************************************************************
V = (V2-V1)/2;
gamma = (gamma2-gamma1)/2;
chi = (chi2-chi1)/2;
h = (h2-h1)/2;
phi = (phi2-phi1)/2;
theta = (theta2-theta1)/2;
mf = (mf2-mf1)/2;

%*********************************************************************
% Sub-Functions
%*********************************************************************

% Alfa Equation
% syms alfav
% dCL1 = 0.8671 + 7.1112*alfav;
% dCD1 = 0.1031 + 4.2706*alfav;
% % Alpha Equation f(alpha)
% f = - psiV*(dCD1) + psigamma*(dCL1/V);
[alfa]=5/57.3;

% Aerodynamic File
H = 10351.8-(0.0368512).*h-(1.02368e-5).*h.^2+(2.63363e-10).*h.^3;
rho = rhos.*exp(-h./H);
g = gs*(rearth/(h + rearth))^2;
beta = 0;

CL = -0.0057 + 0.8671*alfa + 3.5556*alfa^2; %alfa in rad
CD = -0.0101 + 0.1031*alfa + 2.1353*alfa^2;
CQ =  0.0000 - 0.1038*beta + 0.0299*beta^2;

F = 0.5*rho*V^2*Sref*[CD CQ CL];
D = F(1); Y = F(2); L = F(3);
F_v = rho*V*Sref*[CD CQ CL];
Dv = F_v(1); Yv = F_v(2); Lv = F_v(3);

% Thrust File
M = V/a;
I  = 3574.2 + 437.7*M - 401.4*M^2 + 28.8*M^3; % eta = 1.5
dI1  = 437.7 - 802.8*M + 86.4*M^2;
CT = 0.0005*M^4 - 0.0221*M^3 + 0.2763*M^2 - 0.8517*M + 1.0619;
CTd = 0.002*M^3 + 0.0663*M^2 + 0.5526*M - 0.8517;

Thrust  = 0.029*rho*g*a*M*I*Ac*CT;
Thrustd = 0.029*rho*g*a*M*dI1*Ac*CT + 0.029*rho*a*M*g*I*Ac*CTd +... 
          0.029*rho*g*a*I*Ac*CT;

% ********************************************************
% Estimation of Position
% ********************************************************

E4 = h2     - h1 - dt*(V*sin(gamma));
E5 = phi2   - phi1 - dt*(V*cos(gamma)*cos(chi)/(h + rearth));
E6 = theta2 - theta1 - dt*(V*cos(gamma)*sin(chi)/((h + rearth)*cos(phi)));

% ***************************************************************
% Mass Consuming Equation 
% ***************************************************************

E7 = mf2 - mf1 - dt*(-Thrust/(g*I));
m = m0 - mf; 

% ***************************************************************
% Set up of Right Hand Side Eauations
% ***************************************************************
R = h + rearth;
E1 = V2 - V1 - dt*(((Thrust-D)/m)-g*sin(gamma)+omegae^2*R*cos(phi)^2*...
             (sin(gamma)-cos(gamma)*tan(phi)*sin(chi)));
E2 = gamma2 - gamma1 - dt*((1/V)*((L*cos(mu)-Y*sin(mu))/m - g*cos(gamma)+...
        (V^2*cos(gamma))/R+2*omegae*V*cos(phi)*cos(chi)+omegae^2*R*...
         cos(phi)^2*(cos(gamma)+sin(gamma)*tan(phi)*sin(chi))));
E3 = chi2 - chi1 - dt*((1/V)*((L*sin(mu)+Y*cos(mu))/(m*cos(gamma))-(V^2/R)...
    *cos(gamma)*cos(chi)*tan(phi)+2*omegae*V*(tan(gamma)*cos(phi)*...
    sin(chi)-sin(phi))-(omegae^2*R/cos(gamma))*sin(phi)*cos(phi)*cos(chi)));

% Right-Hand Side
En = -[E1;E2;E3;E4;E5;E6;E7];
E0 = -[h2 - 30000;
       phi2 - 33/57.3;
       theta2 - 55/57.3;
       mf2 - 80000];
   
EM = -[h2 - 30000;
       phi2 - 45/57.3;
       theta2 + 75/57.3];

