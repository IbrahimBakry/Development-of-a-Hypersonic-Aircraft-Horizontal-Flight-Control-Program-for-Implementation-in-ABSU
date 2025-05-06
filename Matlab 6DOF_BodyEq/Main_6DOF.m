% ********************************************************
% Main Program 
% Integrating of Translational and Rotational Governing 
% Equations of Flight Mechanics; Six Degrees of Freedom
% Body Fixed Coordinates.
% ********************************************************
clc; clear;

tstart=clock;

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,...
    CLX,CDX,b1X,CyX,CnX]=constants;

% ********************************************************
% Time Parameter
% ********************************************************

arc   = pi/180.;
t0    = 0.;
tend1 = 10.;
tspan = [t0 tend1];

% ********************************************************
% Initial Conditions
% ********************************************************

[un,vn,wn,phisn,thetasn,psisn,hn,phin,thetan]=initiala;

[omb,ome]=initialb(un,vn,wn,phisn,thetasn,psisn,hn,phin,thetan,omegae);

pn = omb(1) + ome(1);
qn = omb(2) + ome(2);
rn = omb(3) + ome(3);

psiu = 0.01;      psiw = 0.01; psiqt = 0.01;
psithetae = 0.01; psih = 0.01; psiphi = 0.01;
psim = 0.01;

% ********************************************************
% Integration Loop
% ********************************************************
    
y0(1)=  un;
y0(2)=  vn;
y0(3)=  wn;
y0(4)=  pn;
y0(5)=  qn;
y0(6)=  rn;
y0(7)=  phisn;
y0(8)=  thetasn;
y0(9)=  psisn;
y0(10)= hn;
y0(11)= phin;
y0(12)= thetan;
y0(13)= psiu;
y0(14)= psiw;
y0(15)= psiqt;
y0(16)= psithetae;
y0(17)= psih;
y0(18)= psiphi;
y0(19)= psim;


options = odeset('Maxstep',1);

[t,y]=ode23(@RHS,tspan,y0,options);
% [t,y]=RK2(tspan,y0);

un      = y(:,1);
vn      = y(:,2);
wn      = y(:,3);
pn      = y(:,4); 
qn      = y(:,5);
rn      = y(:,6);
phisn   = y(:,7);
thetasn = y(:,8);
psisn   = y(:,9);
hn      = y(:,10);
phin    = y(:,11);
thetan  = y(:,12);
psiu      = y(:,13);
psiw      = y(:,14);
psiqt     = y(:,15);
psithetae = y(:,16);
psih      = y(:,17);
psiphi    = y(:,18);
psim      = y(:,19);

% ********************************************************
% Evaluation of Data and Plot Preparation
% ********************************************************

elem = size(un);
iend = elem(1);

nn = 1;
j = 0;
for i = 1:4:iend
  j = j + 1;

V1(1)  = un(i); V1(2)  = vn(i); V1(3)  = wn(i);
V1g    = Mg'*V1';

gamma(j)   =-asin(V1g(3)/Vabs);
chi(j)     = atan(V1g(2)/V1g(1)); 

if V1g(1) < 0, chi(j) = pi + atan(V1g(2)/V1g(1)); end  

[Maero,Vabs,Mg,alfa(1),beta]...
=matrices(un(i),vn(i),wn(i),phisn(i),thetasn(i),psisn(i));

rgn(j) = hn(i) + rearth;

[A,M,MG1,CL,CD,Cm,rho,g,CA,ma]...
=aerodynamics(hn(i),Maero,Vabs,rgn(j),Mg,alfa(j),beta);

alfa(j+1) = double(alpha_equation(Vabs,rho,g,gamma(j),rgn(j),phin(j),...
    psiu(j),psiw(j),psiqt(j),alfa(j),Maero,beta));
% alfa limits
arc = pi/180; 
if alfa(j+1) <= -3*arc, alfa(j+1) = -3*arc; end
if alfa(j+1) >= 21*arc, alfa(j+1) = 21*arc; end


alfagr(j)  = alfa(j+1)/arc;
betagr(j)  = beta/arc;
thetsgr(j) = thetasn(i)/arc; % Euler
psisgr(j)  = psisn(i)/arc;   % Euler
thetgr(j)  = thetan(i)/arc;
phigr(j)   = phin(i)/arc;

gammagr(j) = gamma(j)/arc;
chigr(j)   = chi(j)/arc;
Vmag(j)    = Vabs;
h(j)       = hn(i);

result(nn,1) = t(i);
result(nn,2) = Vmag(j);
result(nn,3) = gammagr(j);
result(nn,4) = alfagr(j);
result(nn,5) = thetsgr(j);
result(nn,6) = psisgr(j);
result(nn,7) = thetgr(j);
result(nn,8) = phigr(j);
result(nn,9) = h(j);
result(nn,10) = betagr(j);
result(nn,11) = chigr(j);
% result(nn,12) = phin(i);
% result(nn,13) = thetan(i);

% result1(nn,1) = t(i);
% result1(nn,2) = CL(i);
% result1(nn,3) = CD(i);
% result1(nn,4) = Cm(i);

nn = nn+1;

end

% ********************************************************
% Generation of Files for Plotting the Results
% ********************************************************
subplot(421)
plot(result(:,1),result(:,2)); title('Vmag')  %Vmag
subplot(422)
plot(result(:,1),result(:,9)); title('h')%h
subplot(423)
plot(result(:,1),result(:,4)); title('alphagrad') %alphagrad
subplot(424)
plot(result(:,1),result(:,10));title('betagrad')%betagrad
subplot(425)
plot(result(:,1),result(:,3)); title('gamma')%gamma
subplot(426)
plot(result(:,1),result(:,11));title('chi')%chi
subplot(427)
plot(result(:,1),result(:,8));title('phi')%phi
subplot(428)
plot(result(:,1),result(:,7));title('theta')%theta

% ********************************************************
% End of Main Program
% ********************************************************