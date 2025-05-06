function dydt= RHS(t,y)
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =
% Sub Program
% Providing the Right H and Sides of the Differential
% Equation dy/dt= by Gauss Elimination of A y(t) = B .
% Actual State: E quations with Earth Rotation Terms
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =

[m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,a1X,CmX,CLX,CDX,b1 X,CyX,CnX]=constants;
[T] = inertia;
u1 = y(1);
v1 = y(2);
w1 = y(3);
p1 = y(4);
q1 = y(5);
r1 = y(6);
phis1 = y(7);
thetas1 = y(8);
psis1 = y(9);
h1 = y(10);
phi1 = y(11);
theta1 = y(12);
[Maero,Vabs,Mg,alfa,beta]= matrices(u1,v1,w1,phis1,thetas1,psis1);
rg1 = h1 + rearth;
[A,M,MG1,CL,CD,Cm]= aerodynamics(h1,Maero,Vabs,rg1,Mg,alfa,beta);

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
% E sta blishment o f E limina tio n Ma tr ix AA
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
AA = zeros(6,6);
AA(1,1)= m;
AA(2,2)= m;
AA(3,3)= m;
AA(4,4)= T(1,1);
AA(4,6)= T(1,3);
AA(5,5)= T(2,2);
AA(6,4)= T(3,1);
AA(6,6)= T(3,3);
om1(1) = cos(phi1) * omegae;
om1(2) = 0.;
om1(3) = - sin(phi1) * omegae ;
omb = Mg * om1';
pb = omb(1);
qb = omb(2);
rb = omb(3);
ommat= zeros(3,3);
ommat(1,2) = -rb;
ommat(1,3) = qb;
ommat(2,1) = rb;
ommat(2,3) = -pb;
ommat(3,1) = -qb;
ommat(3,2) = pb;
ommat2 = ommat * ommat;
rg1v(1) = 0.;
rg1v(2) = 0.;
rg1v(3) = rg1;
rg1vb = Mg * rg1v';
centripetal = -ommat2 * rg1vb;
V1(1) = u1;
V1(2) = v1;
V1(3) = w1;
V1g = Mg'*V1';
gamma = -asin(V1g(3)/Vabs);
chi = atan(V1g(2)/V1g(1));

  if V1g (1) < 0.
    chi = pi + atan(V1g(2)/V1g(1));
  end
h1p = Vabs*sin(gamma);
phi1p = Vabs*cos(gamma)*cos(chi)/rg1;
theta1p = Vabs*cos(gamma)*sin(chi)/(rg1*cos(phi1));
om2(1) = theta1p*cos(phi1);
om2(2) =-phi1p;
om2(3) =-theta1p*sin(phi1);
ome = Mg*om2';
pe = ome(1);
qe = ome(2);
re = ome(3);

% ***************************************************************
% Set up of Right Hand Side (The Inhomogeneous part of the ODE)
% ***************************************************************
BB(1) = A(1) + MG1(1) - m*(q1*w1-r1*v1) - m*(qb*w1-rb*v1) - m*centripetal(1);
BB(2) = A(2) + MG1(2) - m*(r1*u1- p1*w1) - m*(rb*u1-pb*w1) - m*centripetal(2);
BB(3) = A(3) + MG1(3) - m*(p1*v1- q1*u1) - m*(pb*v1-qb*u1) - m*centripetal(3);
BB(4)=M(1)-T(1,3)*p1*q1+(T(2,2)-T(3,3))*q1*r1;
BB(5)=M(2)-T(1,3)*(r12-p12)+ (T(3,3)-T(1,1))*r1*p1;
BB(6)=M(3)-T(1,3)*q1*r1+(T(1,1)-T(2,2))*p1*q1;
X=inv(AA)*BB';
pp1 = p1 - pb - pe;
qq1 = q1 - qb - qe;
rr1 = r1 - rb - re;

% ***************************************************************
% Rate of Change of the Euler Angles Describing the Attitude
% of the Space Vehicle
% ***************************************************************
phis1p = pp1+sin(phis1)*tan(thetas1)*qq1+cos(phis1)*tan(thetas1)*rr1;
thetas1p = cos(phis1)*qq1-sin(phis1)*rr1;
psis1p = sin(phis1)/cos(thetas1)*qq1+cos(phis1)/cos(thetas1)*rr1;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** *
% Completion of the Set of the 12 Ordinary Differential Equations
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** *
dydt= X;
dydt(7) = phis1p;
dydt(8) = thetas1p;
dydt(9) = psis1p;
dydt(10) = h1p;
dydt(11) = phi1p;
dydt(12) = theta1p;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **
% End of Sub Program RHS
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **