clc; clear all

a = 0; b = 1;
N = 50; dx = (b-a)/N; x = (a:dx:b)';
y1(1) = 1; y2(1)=0;
for k = 1:N+1
if k >= 2
b(2*k+1:2*k+2,1) = 0;
b(1,1)=1;
b(2*k+2,1)=3;
b(2*k:2*k+1,1) = [dx*(4*(y1(k))*(y2(k))-(y1(k))^2);
                  dx*((y1(k))^2+(y2(k))^2)];
A(2*k+1:2*k+2,2*k+1:2*k+2) = 0;
A(1,1) = 1.0; 
A(2*k:2*k+1,2*k-1:2*k+2) = [-1  0 1 0;
                             0 -1 0 1];
A(2*k+2,2*k+2) = 1.0;
% solving the matrix problem looks like:
y = A \ b; % solve A Y = b

% y1(k)= y1(k-1) + dx*y(2*k-1); y2(k)= y2(k-1) + dx*y(2*k);
% y1(k+1)= y1(k) + dx*y(2*k+1); y2(k+1)= y2(k) + dx*y(2*k+2);
% y1(k) = y(2*k-1); y2(k) = y(2*k);
y1(k+1) = y(2*k+1); y2(k+1) = y(2*k+2);
else
b = zeros(2*k+2,1);
b(1,1)=1;  
b(2*k+2,1)=3;
b(2*k:2*k+1,1) = [dx*(4*y1(k)*y2(k)-y1(k)^2);
                  dx*(y1(k)^2+y2(k)^2)];
A = zeros(2*k+2,2*k+2);
A(1,1) = 1.0; 
A(2*k:2*k+1,2*k-1:2*k+2) = [-1  0 1 0;
                             0 -1 0 1];
A(2*k+2,2*k+2) = 1.0;
% solving the matrix problem looks like:
y = A \ b; % solve A Y = b
y1(k+1)= y(2*k+1); y2(k+1)=y(2*k+2);
end
end

% also get exact soln on fine grid:
subplot(211);plot(x,y1(1:end-1),'-o','markersize',4);hold on;grid on
subplot(212);plot(x,y2(1:end-1),'-o','markersize',4);hold on;grid on
