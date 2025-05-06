% Find u over t in [0; t ] to maximize
% J = mf 
clc; clear
% Problem setup
toms t t_f

% Scaled time
p1 = tomPhase('p1', t, 0, t_f, 300);
setPhase(p1);
tomStates V gamma chi h x y mf 
tomControls alpha throttle

% Constants
m0     = 136077.7;
omegae = 0.00007292;
rearth = 6.38E6;
alt0   = 30000;         altf   = 30000;
vel0   = 1830;          velf   = 1830;
tMin   = 100;           tMax   = 4000;
tGuess = 2000;
arc    = pi/180;
rho = 0.0184;
grav = 9.7151;
beta = 0;
mu = 0;
% Initial guess
x0 = {
t_f == 100
icollocate({
h == alt0-(alt0-altf)*t/t_f
x == 0
y == 0
V == vel0-(vel0-velf)*t/t_f
gamma == -1*arc-2*arc*t/t_f
chi == 1*arc-1*arc*t/t_f
mf == 9
})
collocate({
alpha == 0.03
throttle == 0.5
})
};

% Initial and Final constraints
cbnd = {
initial({
h == alt0
x == 0
y == 0
V == 1830
gamma == 1*arc
chi == 1*arc
mf == 9
alpha == 0.03
throttle == 0.5
})
final({
h == altf
V == velf
x == 30000
y == 50000
})};

% Box constraints
cbox = {
500 <= t_f <= 5000
29995 <= icollocate(h) <= 30005
0 <= icollocate(x) <= 30000
0 <= icollocate(y) <= 50000
1825 <= icollocate(V) <= 1835
-45*arc <= icollocate(gamma) <= 45*arc
-89*arc <= icollocate(chi) <= 89*arc
9.071 <= collocate(mf) <= 163.29
-15*arc <= collocate(alpha) <= 15*arc
0 <= collocate(throttle) <= 1
};
% Additional Algebric Equations
D = 0.1994.*V.^2.*(-0.0101 + 0.1031.*alpha + 2.1353.*alpha.^2);
L = 0.1994.*V.^2.*(-0.0057 + 0.8671.*alpha + 3.5556.*alpha.^2);
Y = 0.1994.*V.^2.*(+0.0000 - 0.1038.*beta + 0.0299.*beta.^2);
Thrust  = grav.*(9.071+154.219.*throttle)*(3574.2+1.4351.*V-0.0043.*V.^2+28.8*(V/305)^3);

% ODEs and path constraints
ceq = collocate({
dot(V)      == ((Thrust-D)./(m0-mf))-grav.*sin(gamma)
dot(gamma)  == (1/V).*((L.*cos(mu)-Y.*sin(mu))./(m0-mf) - grav.*cos(gamma))
dot(chi)    == (1/V).*((L.*sin(mu)+Y.*cos(mu))./((m0-mf).*cos(gamma)))
dot(h)      == V.*sin(gamma)
dot(x)      == V.*cos(gamma).*cos(chi)
dot(y)      == V.*cos(gamma).*sin(chi)
dot(mf)     == -(9.071+(163.29-9.071).*throttle)
});

% Objective
objective = -final(mf);

% Solve the problem
options = struct;
options.name = 'Level Flight';
options.Prob.SOL.optPar(30) = 100000;
options.scale = 'auto';
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);
% 
% for i = 1:length(mf)
%     Mf(i) = m0 - mf(i);
% end

% Plot result
close all
figure
subplot(311); ezplot(h); legend('alt'); title('Altitude');
subplot(312); ezplot(V); legend('vel'); title('Velocity');
subplot(313); ezplot(mf); legend('mf'); title('Mf');
figure
subplot(411); ezplot(x); legend('x'); title('X');
subplot(412); ezplot(y); legend('y'); title('Y');
subplot(413); ezplot(gamma); legend('gamma'); title('Gamma');
subplot(414); ezplot(chi); legend('chi'); title('Chi');
figure
subplot(311); ezplot(alpha); legend('alpha'); title('Alpha');
% subplot(312); ezplot(mu); legend('mu');title('Beta');
subplot(313); ezplot(throttle); legend('throttle');title('Throttle');
