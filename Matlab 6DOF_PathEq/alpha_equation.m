function [alfo]=alpha_equation(V,psiV,psigamma,psiqt)
% Solving Alpha Equation  

[m,Sref,Lref] = constants;
[T] = inertia; Iy = T(2,2);

syms alfav
dCL1 = 0.8671 + 7.1112*alfav;
dCD1 = 0.1031 + 4.2706*alfav;
dCm1 = 0.1384 - 1.5258*alfav;

dCL2 =  7.1112;
dCD2 =  4.2706;
dCm2 = -1.5258;

% Alpha Equation f(alpha)
f = - psiV*(dCD1) + psigamma*(dCL1/V) + psiqt*(dCm1*(m*Lref/Iy));

% Equation f_dot(alpha)
f_dot = - psiV*(dCD2) + psigamma*(dCL2/V) + psiqt*(dCm2*(m*Lref/Iy));

N = 1000; % ? of Iterations
err = 0.001; x0 = 0.1;
alfo = NRaphson(f,f_dot,N,err,x0);
arc = pi/180;
if alfo <= (-03)*arc, alfo = -03*arc; end
if alfo >= (+21)*arc, alfo = +21*arc; end