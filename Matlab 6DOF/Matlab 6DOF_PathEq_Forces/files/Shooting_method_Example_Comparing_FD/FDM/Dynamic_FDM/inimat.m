function [S_k1,S_k,S0,SM] = inimat(k,dt,y)
% deriving dervitevs
syms V1 V2 gamma1 gamma2 chi1 chi2 h2 h1 phi2 phi1 theta2 theta1 mf2 mf1
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


% Main Matrices


S_k1 = [diff(E1,V1) diff(E1,gamma1) diff(E1,chi1) diff(E1,h1) diff(E1,phi1) diff(E1,theta1) diff(E1,mf1);
        diff(E2,V1) diff(E2,gamma1) diff(E2,chi1) diff(E2,h1) diff(E2,phi1) diff(E2,theta1) diff(E2,mf1);
        diff(E3,V1) diff(E3,gamma1) diff(E3,chi1) diff(E3,h1) diff(E3,phi1) diff(E3,theta1) diff(E3,mf1);
        diff(E4,V1) diff(E4,gamma1) diff(E4,chi1) diff(E4,h1) diff(E4,phi1) diff(E4,theta1) diff(E4,mf1);
        diff(E5,V1) diff(E5,gamma1) diff(E5,chi1) diff(E5,h1) diff(E5,phi1) diff(E5,theta1) diff(E5,mf1);
        diff(E6,V1) diff(E6,gamma1) diff(E6,chi1) diff(E6,h1) diff(E6,phi1) diff(E6,theta1) diff(E6,mf1);
        diff(E7,V1) diff(E7,gamma1) diff(E7,chi1) diff(E7,h1) diff(E7,phi1) diff(E7,theta1) diff(E7,mf1)];
S_k1 =  double(subs(S_k1,[V1 gamma1 chi1 h1 phi1 theta1 mf1 V2 gamma2 chi2 h2 phi2 theta2 mf2]...
                        ,[0 0 0 0 0 0 0 y(k,1) y(k,2) y(k,3) y(k,4) y(k,5) y(k,6) y(k,7)]));

S_k = [diff(E1,V2) diff(E1,gamma2) diff(E1,chi2) diff(E1,h2) diff(E1,phi2) diff(E1,theta2) diff(E1,mf2);
       diff(E2,V2) diff(E2,gamma2) diff(E2,chi2) diff(E2,h2) diff(E2,phi2) diff(E2,theta2) diff(E2,mf2);
       diff(E3,V2) diff(E3,gamma2) diff(E3,chi2) diff(E3,h2) diff(E3,phi2) diff(E3,theta2) diff(E3,mf2);
       diff(E4,V2) diff(E4,gamma2) diff(E4,chi2) diff(E4,h2) diff(E4,phi2) diff(E4,theta2) diff(E4,mf2);
       diff(E5,V2) diff(E5,gamma2) diff(E5,chi2) diff(E5,h2) diff(E5,phi2) diff(E5,theta2) diff(E5,mf2);
       diff(E6,V2) diff(E6,gamma2) diff(E6,chi2) diff(E6,h2) diff(E6,phi2) diff(E6,theta2) diff(E6,mf2);
       diff(E7,V2) diff(E7,gamma2) diff(E7,chi2) diff(E7,h2) diff(E7,phi2) diff(E7,theta2) diff(E7,mf2)];
S_k = double(subs(S_k,[V1 gamma1 chi1 h1 phi1 theta1 mf1 V2 gamma2 chi2 h2 phi2 theta2 mf2]...
                     ,[0 0 0 0 0 0 0 y(k,1) y(k,2) y(k,3) y(k,4) y(k,5) y(k,6) y(k,7)]));

S0 = [0 0 0 1 0 0 0;  % h0,phi0,theta0,mf0
      0 0 0 0 1 0 0;
      0 0 0 0 0 1 0;
      0 0 0 0 0 0 1];
SM = [0 0 0 1 0 0 0; % hf,phif,thetaf
      0 0 0 0 1 0 0;
      0 0 0 0 0 1 0];