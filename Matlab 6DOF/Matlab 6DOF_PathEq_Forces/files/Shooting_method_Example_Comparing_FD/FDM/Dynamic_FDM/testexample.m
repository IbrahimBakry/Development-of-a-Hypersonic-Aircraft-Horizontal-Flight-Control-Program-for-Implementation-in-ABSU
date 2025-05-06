clc; clear all

a = 0; b = 20;
N = 50; dx = (b-a)/N; x = (a:dx:b)';
y(1) = 1;
n = 1; cf = 1;
for k = 1:N+1

b(n*k:n*k+cf,1) = 0;
b(n*k,1) = [dx*2*y(k)^2];
b(n*k+cf,1) = 2;

A(n*k+1:n*k+1,n*k+1:n*k+1) = 0;
A(n*k,n*k:n*k+1) = [-1  1];
A(n*k+1,n*k+1) = 1.0;

% solving the matrix problem looks like:
S = A \ b; % solve A Y = b

% y(k) = y(k-1) + dx*S(n*k);
% y(k+1) = y(k) + dx*S(n*k+1);
% y(k)   = S(n*k);
y(k+1) = S(n*k+1);

end

% also get exact soln on fine grid:
subplot(211);plot(x,y(1:end-1),'-o','markersize',4);hold on;grid on
