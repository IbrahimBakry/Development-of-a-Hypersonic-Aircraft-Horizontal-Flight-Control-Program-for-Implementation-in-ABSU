clc; clear all;
% bvp4c Example
% solinit = bvpinit(linspace(0,1,4),[1 0 0 0]);
% sol = bvp4c(@twoode,@twobc,solinit);
% x = linspace(0,1);
% y = deval(sol,x);
% plot(x,y(1,:))
% subplot(211);plot(x,y(1,:)); grid on; title('y1')
% subplot(212);plot(x,y(2,:)); grid on; title('y2')

% bvp4c Dynamica
% solinit = bvpinit(linspace(0,100,100),[1830,0.5,0.5,30000,33/57.3,55/57.3,80000]);
% sol = bvp4c(@dynamic2,@dtwobc,solinit);
% x = linspace(0,100);
% y = deval(sol,x);
% subplot(421);plot(x,y(1,:)); grid on; title('V')
% subplot(422);plot(x,y(2,:).*57.3); grid on; title('gamma')
% subplot(423);plot(x,y(3,:).*57.3); grid on; title('chi')
% subplot(424);plot(x,y(4,:)); grid on; title('h')
% subplot(425);plot(x,y(5,:).*57.3); grid on; title('phi')
% subplot(426);plot(x,y(6,:).*57.3); grid on; title('theta')
% subplot(427);plot(x,y(7,:)); grid on; title('mf')

solinit = bvpinit(linspace(0,100,150),[1 0]);
sol = bvp4c(@twoode,@twobc,solinit);
x = linspace(0,100);
y = deval(sol,x);
subplot(211);plot(x,y(1,:));hold on
subplot(212);plot(x,y(2,:));hold on


