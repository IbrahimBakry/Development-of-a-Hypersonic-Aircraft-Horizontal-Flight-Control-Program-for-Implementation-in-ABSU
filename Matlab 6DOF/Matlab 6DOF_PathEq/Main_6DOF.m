% ********************************************************
% Main Program 
% Integrating of Translational and Rotational Governing 
% Equations of Flight Mechanics; Six Degrees of Freedom
% Body Fixed Coordinates.
% ********************************************************
clc; clear;

tstart=clock;

[m,Sref,Lref,bref,Ac,a,rhos,gs,bet,omegae,rearth,V00]=constants;

% ********************************************************
% Time Parameter
% ********************************************************

arc   = pi/180.;
t0    = 0.;
tend1 = 20.;
tspan = [t0 tend1];

% ********************************************************
% Initial Conditions
% ********************************************************

[Vn,gamman,chin,phisn,thetasn,psisn,hn,phin,thetan,mfn] = initiala;

[omb,ome] = initialb(phisn,thetasn,psisn,hn,phin,thetan,omegae);

pn = omb(1) + ome(1);
qn = omb(2) + ome(2);
rn = omb(3) + ome(3);

psiV = 0.0001;   psigamma = 0.0001;   psichi = 0.0001; 
psiqt = 0.0001;  psithetae = 0.0001;  psih = 0.0001; 
psiphi = 0.0001; psitheta = 0.0001;   psim = 0.0001;       

% psiV2     = 0.01; psigamma2 = 0.01; psichi2   = 0.01;
% psipt2    = 0.01; psirt2    = 0.01; psiphie2  = 0.01;
% psipsie2  = 0.01; psiphi2   = 0.01; psitheta2 = 0.01;
% psim2     = 0.01;

% ********************************************************
% Integration Loop
% ********************************************************
    
y0(1)=  Vn;         y0(2)=  gamman;   y0(3)=  chin;
y0(4)=  pn;         y0(5)=  qn;       y0(6)=  rn;
y0(7)=  phisn;      y0(8)=  thetasn;  y0(9)=  psisn;
y0(10)= hn;         y0(11)= phin;     y0(12)= thetan;
y0(13)= mfn;        y0(14)= psiV;     y0(15)= psigamma;
y0(16)= psichi;     y0(17)= psiqt;    y0(18)= psithetae;
y0(19)= psih;       y0(20)= psiphi;   y0(21)= psitheta;
y0(22)= psim;      % y0(23)= psiV2;    y0(24)= psigamma2;
% y0(25)= psichi2;    y0(26)= psipt2;   y0(27)= psirt2;
% y0(28)= psiphie2;   y0(29)= psipsie2; y0(30)= psiphi2;
% y0(31)= psitheta2;  y0(32)= psim2;


options = odeset('Maxstep',1);

[t,y]=ode23(@RHS,tspan,y0,options);
% [t,y]=RK2(tspan,y0);

Vn      = y(:,1);
gamman  = y(:,2);
chin    = y(:,3);
pn      = y(:,4); 
qn      = y(:,5);
rn      = y(:,6);
phisn   = y(:,7);
thetasn = y(:,8);
psisn   = y(:,9);
hn      = y(:,10);
phin    = y(:,11);
thetan  = y(:,12);
mfn      = y(:,13);
psiV      = y(:,14);
psigamma  = y(:,15);
psichi    = y(:,16);
psiqt     = y(:,17);
psithetae = y(:,18);
psih      = y(:,19);
psiphi    = y(:,20);
psitheta  = y(:,21);
psim      = y(:,22);


% % ********************************************************
% % Generation of Files for Plotting the Results
% % ********************************************************
subplot(421)
plot(t,Vn); title('Vmag')  %Vmag
subplot(422)
plot(t,hn); title('h')%h
subplot(423)
plot(t,mfn); title('mf') %mf
% subplot(424)
% plot(result(:,1),result(:,10));title('betagrad')%betagrad
subplot(425)
plot(t,gamman.*57.3); title('gamma')%gamma
subplot(426)
plot(t,chin.*57.3);title('chi')%chi
subplot(427)
plot(t,phin.*57.3);title('phi')%phi
subplot(428)
plot(t,thetan.*57.3);title('theta')%theta
% 
% % ********************************************************
% % End of Main Program
% % ********************************************************