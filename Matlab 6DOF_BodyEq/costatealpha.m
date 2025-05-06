function dydt=costatealpha(t,y)
%%% Equations of Costate of alpha

% rho, u, w, V, Cx, Cz, Cm, qt, R, 
% g, thetae, phi, pe, re, gamma, Tx

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,...
          CLX,CDX,b1X,CyX,CnX]=constants;
[T] = inertia; Iy = T(2,2);

psiu = y(1);
psiw = y(2);
psiqt = y(3);
psithetae = y(4);
psih =  y(5);
psiphi = y(6);
psim = y(7);

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
    
psiphid = -psiu(R*omegae^2*(2*cos(thetae)*cos(phi)^2+2*cos(phi)...
    *sin(thetae)*sin(phi)-cos(thetae))) - psiw*(R*omegae^2*(2*sin(thetae...
    )*cos(phi)^2-2*cos(phi)*cos(thetae)*sin(phi)-sin(thetae)));

psimd = psiu*((0.5*rho*V^2*Sref*Cx+Tx)/m^2) + ...
    psiw*((0.5*rho*V^2*Sref*Cz)/m^2);



dydt(1) = psiud;
dydt(2) = psiwd;
dydt(3) = psiqtd;
dydt(4) = psithetaed;
dydt(5) = psihd;
dydt(6) = psiphid;
dydt(7) = psimd;







