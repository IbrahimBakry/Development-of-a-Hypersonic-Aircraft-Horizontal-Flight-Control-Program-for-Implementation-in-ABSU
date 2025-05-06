function dyt = dynamic2(t,y)
%*********************************************************************
% Constants
%*********************************************************************
m0      = 136077.7;% Takeoff weight;
omegae = 0.00007292;
rearth = 6.38E6;

%*********************************************************************
% Assigning Variables
%*********************************************************************
V1      = y(1);    gamma1  = y(2);   chi1    = y(3);
h1      = y(4);    phi1    = y(5);   theta1  = y(6);
mf      = y(7);   

[alfa]=5/57.3;
[D,Y,L,Dv,Yv,Lv,g,rho] = aero([V1 h1 alfa]);
[Thrust,Thrustd,I] = thrust([V1 g rho]);

F = [D Y L];      F_v = [Dv Yv Lv]; 

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
dydt = zeros(7,1);
dydt(1,1) = x1(1);
dydt(2,1) = x1(2);
dydt(3,1) = x1(3);

dydt(4,1) = h1p;
dydt(5,1) = phi1p;
dydt(6,1) = theta1p;  
dydt(7,1) = mfd;

dyt = double(dydt);

% ***************************************************************
% End of Sub Program RHS
% ***************************************************************
