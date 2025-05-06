function [Thrust, Thrustd, I] = thrust(u)

Ac     = 27.87;
a      = 305;

V1 = u(1);
g = u(2);
rho = u(3);

Mach = V1/a;
[I,dI1,dI2,p,pd]=Engine_CT(Mach);

CTP  = polyval(p,Mach);
CTPd = polyval(pd,Mach);

Thrust  = 0.029*rho*g*V1*I*Ac*CTP;
Thrustd = 0.029*rho*g*dI1*Ac*CTP*V1 + 0.029*rho*g*I*Ac*CTPd*V1 +... 
          0.029*rho*g*(1/a)*I*Ac*CTP;


