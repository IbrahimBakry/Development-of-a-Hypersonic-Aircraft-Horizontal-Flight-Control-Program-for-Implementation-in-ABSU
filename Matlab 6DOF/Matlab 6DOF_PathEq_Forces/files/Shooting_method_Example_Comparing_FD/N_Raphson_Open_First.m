clc
clear all
close all

n=1000;
alpha=0.1;
yb=3;
Error = 1;

while Error >= 1e-3
[xx,zz]=RK4shoot(n,alpha);
 alpha = alpha - (zz(n,2) - yb)/zz(n,3)

 Error = abs(zz(n,2) - yb)
end

disp({'alpha = ',alpha})
subplot(211);plot(xx,zz(:,1)); grid on; title('y1')
subplot(212);plot(xx,zz(:,2)); grid on; title('y2')


 