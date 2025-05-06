% =========================================
% Main Program
% Integrating of Translational and Rotational Governing
% Equations of Flight Mechanics; Six Degree of Freedom
% Body Fixed Coordinates.
% =========================================
tstart=clock;
[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,a1X,CmX,CLX,CDX,b1X,CyX,CnX]=constants;

% ********************************************************
% Time Parameter
% ********************************************************
arc = pi/180.;
t0 = 0.;

tend1 = 989.;
tspan = [t0 tend1];

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Initial Conditions
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
[un,vn,wn,phisn,thetasn,psisn,hn,phin,thetan]= initiala ;
[omb,ome]= initialb(un,vn,wn,phisn,thetasn,psisn,hn,phin,thetan,omegae);
pn = omb(1) + ome(1);
qn = omb(2) + ome(2);
rn = omb(3) + ome(3);

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Integration Loop
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
y0(1)= un;
y0(2)= vn;
y0(3)= wn;
y0(4)= pn;
y0(5)= qn;
y0(6)= rn;
y0(7)= phisn;
y0(8)= thetasn;
y0(9)= psisn;
y0(10)= hn;
y0(11)= phin;
y0(12)= thetan;
options = odeset('Maxstep',1);
[t,y]= ode23(@RHS,tspan,y0,options);
un = y(:,1);
vn = y(:,2);
wn = y(:,3);
pn = y(:,4);
qn = y(:,5);
rn = y(:,6);
phisn = y(:,7);
thetasn = y(:,8);
psisn = y(:,9);
hn = y(:,10);
phin = y(:,11);
thetan = y(:,12);

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Evaluation of Data and Plot Preparation
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
elem = size(un);
iend = elem(1);
nn = 1;
for i = 1:4:iend
[Maero ,Vabs,Mg ,alfa ,beta ]= matrices(un(i),vn(i),wn(i),phisn(i),thetasn(i),psisn(i));
rgn(i) = hn(i) + rearth;
[A,M,MG1,CL(i),CD(i),Cm(i)]= aerodynamics(hn(i),Maero,Vabs,rgn(i),Mg ,alfa ,beta);
V1(1) = un(i);
V1(2) = vn(i);
V1(3) = wn(i);
V1g = Mg'*V1';
gamma(i) = -asin(V1g(3)/Vabs);
chi(i) = atan(V1g(2)/V1g(1));
  if V1g(1) < 0.
    chi(i) = pi + atan(V1g(2)/V1g(1));
  end
alfagr(i) = alfa/arc;
betagr(i) = beta/arc;
thetsgr(i) = thetasn(i)/arc;
thetgr(i) = thetan(i)/arc;
psisgr(i) = psisn(i)/arc;
phigr(i) = phin(i)/arc;

gammagr(i) = gamma(i)/arc;
chigr(i) = chi(i)/arc;
Vmag(i) = Vabs ;
h(i) = hn(i);
result(nn,1) = t(i);
result(nn,2) = Vmag(i);
result(nn,3) = gammagr(i);
result(nn,4) = alfagr(i);
result(nn,5) = thetsgr(i);
result(nn,6) = psisgr(i);
result(nn,7) = thetgr(i);
result(nn,8) = phigr(i);
result(nn,9) = h(i);
result(nn,10)= betagr(i);
result1(nn,1) = t(i);
result1(nn,2) = CL(i);
result1(nn,3) = CD(i);
result1(nn,4) = Cm(i);
nn = nn + 1 ;
end

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% Generation of Files for Plotting the Results
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
fid = fopen('X-38_six_degrees_1.txt','w');
fprintf(fid,'time velocity gamma alfa thetas psis theta phi alt beta \n \n');
fprintf(fid,'%4.0f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %13.4f%10.4f \n',result');
status = fclose(fid);
fid = fopen('X-38_Aero_six_degrees_1.txt','w');
fprintf(fid,'lift drag pitch \n \n');
fprintf(fid,'%4.0f %11.5f %11.5f %11.7f \n',result1');
status = fclose(fid);

% ****************************************************************
% End of Main Program
% ****************************************************************