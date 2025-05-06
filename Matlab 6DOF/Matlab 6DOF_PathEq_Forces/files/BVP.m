clear all;
N=20;   % Number of iterations
beta=2; % Target value
yp=1;   % y’(0) guess
alpha=.1; % Fraction of difference
sim('BVPmodel')
for i=1:N
yp=yp-alpha*(y(end,1)-beta);
sim('BVPmodel')
plot(tout,y,'b')
hold all
end
hold off
axis square
xlabel('x')
ylabel('y')
title('y vs x')