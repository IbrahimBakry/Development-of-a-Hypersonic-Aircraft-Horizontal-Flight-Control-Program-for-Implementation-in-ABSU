function [alfo]=alpha_equation(V,rho,g,gamma,R,phi,psiu,psiw,psiqt,alfo,Maero,beta)
% Solving Alpha Equation  

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,...
    CLX,CDX,b1X,CyX,CnX]=constants;
[T] = inertia; Iy = T(2,2);
arc = pi/180;  qs = 0.5*rho*V^2*Sref;

syms alfav

dCL1 =  0.0150 + 0.0254*alfav + 0.0588*alfav^2 + 0.0594*alfav^3;
dCD1 =  0.0008 - 0.0248*alfav + 0.3720*alfav^2 - 0.6208*alfav^3;
dCm1 = -0.0097 + 0.0672*alfav - 0.4293*alfav^2;
dCy1 = -0.1053 + 0.0530*beta + 0.3024*beta^2;

dCL2 =  0.0254 + 0.1176*alfav + 0.1782*alfav^2;
dCD2 =- 0.0248 + 0.7440*alfav - 1.8624*alfav^2;
dCm2 =  0.0672 - 0.8586*alfav;
dCy2 =  0.0530 + 0.6048*beta;

CXYZ1=Maero*[-dCD1 dCy1 -dCL1]'; dCx1=CXYZ1(1);dCy1=CXYZ1(2);dCz1=CXYZ1(3);
CXYZ2=Maero*[-dCD2 dCy2 -dCL2]'; dCx2=CXYZ2(1);dCy2=CXYZ2(2);dCz2=CXYZ2(3);

% Alpha Equation f(alpha)
f = psiu*(dCx1 - ((m/(qs))*(g*cos(gamma+alfav...
    )+R*omegae^2*(cos(phi)^2*cos(gamma+alfav)+0.5*sin(2*phi)*sin(gamma+...
    alfav))))) + psiw*(dCz1 + ((m/(qs))*(-g*sin(gamma+...
    alfav)+R*omegae^2*(-cos(phi)^2*sin(gamma+alfav)+0.5*sin(2*phi)*...
    cos(gamma+alfav))))) + psiqt*(dCm1*(m*Lref/Iy));

% Equation f_dot(alpha)
f_dot = psiu*(dCx2 + ((m/(qs))*...
    (g*sin(gamma+alfav)+R*omegae^2*(cos(phi)^2*sin(gamma+alfav)-0.5*...
    sin(2*phi)*cos(gamma+alfav))))) + psiw*(dCz2...
    + ((m/(qs))*(-g*cos(gamma+alfav)-R*omegae^2*(cos(phi)^2*...
    cos(gamma+alfav)+0.5*sin(2*phi)*sin(gamma+alfav))))) + psiqt*(...
    dCm2*(m*Lref/Iy));

N = 1000; % ? of Iterations

err = 0.001; x0 = double(alfo);
alfo = NRaphson(f,f_dot,N,err,x0);

%%% Old Formation
% F1 = ((m/(qs))*(g*cos(gamma+alpha)+R*omegae^2*(cos(phi)^2*...
%     cos(gamma+alpha)+0.5*sin(2*phi)*sin(gamma+alpha))));
% F2 = ((m/(qs))*(-g*sin(gamma+alpha)+R*omegae^2*(-cos(phi)^2*...
%     sin(gamma+alpha)+0.5*sin(2*phi)*cos(gamma+alpha))));
% F11 = ((m/(qs))*(g*sin(gamma+alpha)+R*omegae^2*(cos(phi)^2*...
%     sin(gamma+alpha)-0.5*sin(2*phi)*cos(gamma+alpha))));
% F22 = ((m/(qs))*(-g*cos(gamma+alpha)-R*omegae^2*(cos(phi)^2*...
%     cos(gamma+alpha)+0.5*sin(2*phi)*sin(gamma+alpha))));
