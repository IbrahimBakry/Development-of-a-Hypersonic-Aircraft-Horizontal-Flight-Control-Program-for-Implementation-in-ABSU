

x1(1) = 0; x2(1)= 0.1; h=0.01;
ii(1) = 0;
for i = 1:0.1:1
    ii(i+1) = i;
x1(i+1) = x1(i) + h*(x1(i) + x2(i)^2);
x2(i+1) = x2(i) + h*(x1(i)*x2(i) + x2(i));
end

plot(ii,x1,'.-',ii,x2,'.-')


% figure
% [t,y] = ode45(@test0,[0 1],[0 0.1]);
% plot(t,y)
