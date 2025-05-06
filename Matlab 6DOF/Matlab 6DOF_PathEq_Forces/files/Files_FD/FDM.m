% Finite Difference Method

clc; clear all

% Initial Value
Initial_conitions;
y(:,1) = [Vtot0 gamma0 chi0 h0 phi0 theta0 mf0 psiV0 psigamma0...
    psichi0 psih0 psiphi0 psitheta0 psim0]';

for i = 2 : 1000

% Getting Parameters from dynamic
V1      = y(1,i-1);     gamma1  = y(2,i-1);   chi1    = y(3,i-1);
h1      = y(4,i-1);     phi1    = y(5,i-1);   theta1  = y(6,i-1);
mf      = y(7,i-1);     psiV    = y(8,i-1);   psigamma= y(9,i-1);    
psichi  = y(10,i-1);    psih    = y(11,i-1);  psiphi= y(12,i-1);     
psitheta  = y(13,i-1);  psim    = y(14,i-1);
R = h1 + rearth;

[alfa]=alpha_equation(V1,psiV,psigamma);
[D,Y,L,Dv,Yv,Lv,g,rho] = aero([V1 h1 alfa]);
[Thrust,Thrustd,I] = thrust([V1 g rho]);
F = [D Y L];      F_v = [Dv Yv Lv]; 

% CT1 = 0.0005.*(V1/a).^4 - 0.0221.*(V1/a).^3 + 0.2763.*(V1/a).^2 - 0.8517.*(V1/a) + 1.0619;
% CT2 = 0.002.*(V1/a).^3 - 0.0663.*(V1/a).^2 + 0.5526.*(V1/a) - 0.8517;
% mf = mf - (0.029*rho*V1*Ac*CT1)/(0.029*rho*V1*Ac*CT2);
m = m0 - mf;

% Calculating Jacobi
JAK = Jak(F,F_v,g,rho,Thrust,Thrustd,I,m,omegae);
Jako = real(double(subs(JAK)))

% Calculating New Dynamic
dyt = dynamic(y(:,i-1),i-1);

syms V1 gamma1 chi1 h1 phi1 theta1 mf psiV psigamma psichi psiphi psitheta psih psim
DV = [V1 gamma1 chi1 h1 phi1 theta1 mf psiV psigamma psichi psih psiphi psitheta psim]';

E = Jako*DV - Jako*y(:,i-1) + dyt;


[u,it] = sor(Jako,dyt)
% N = 20; err = 0.001; x0 = zeros(14,1);
% x = GSeidelm2(Jako,dyt,x0);
% x = GSeidelm(Jako,dyt,N,err);

% y(:,i) = y(:,i-1) - JAK\(dyt);

end