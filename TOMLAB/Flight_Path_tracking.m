% 36.2 Problem setup
toms t
p = tomPhase('p', t, 0, 100, 60);
setPhase(p);
tomStates x1s x2 x3s
tomControls u1s u2
x1 = x1s*100;
x3 = x3s*10;
u1 = u1s*10e3;
cr2d = pi/180;
% Box constraints
cbox = {
92 <= icollocate(x1) <= 170
-20*cr2d <= icollocate(x2) <= 25*cr2d
-150 <= icollocate(x3) <= 150
60e3 <= collocate(u1) <= 120e3
-150 <= collocate(u2) <= 150};
% Boundary constraints
cbnd = initial({x1 == 153.73; x2 == 0; x3 == 0});
L = 65.3;
D = 3.18;
m = 160e3;
g = 9.81;
c = 6;

% ODEs and path constraints
ceq = collocate({
dot(x1) == (-D/m*x1.^2-g*sin(x2)+u1/m)
dot(x2) == L/m*x1.*(1-c*x2)-g*cos(x2)./x1+L*c/m*u2
dot(x3) == (x1.*sin(x2))});
% Objective
objective = integrate(x3.^2);

% 36.3 Solve the problem
options = struct;
options.name = 'Flight Path Tracking';
solution = ezsolve(objective, {cbox, cbnd, ceq}, [], options);
t = subs(collocate(t),solution);
x1 = subs(collocate(x1),solution);
x2 = subs(collocate(x2),solution);
x3 = subs(collocate(x3),solution);
u1 = subs(collocate(u1),solution);
u2 = subs(collocate(u2),solution);

% 36.4 Plot result
figure(1);
subplot(2,3,1);
plot(t,x3,'-');
title('Alt')
subplot(2,3,2);
plot(t,x1,'-');
title('Vel')
subplot(2,3,3);
plot(t,x2,'-');
title('Gamma')
subplot(2,3,4);
plot(t,u2,'-');
title('Angle')
subplot(2,3,5);
plot(t,u1,'-');
title('Thrust')