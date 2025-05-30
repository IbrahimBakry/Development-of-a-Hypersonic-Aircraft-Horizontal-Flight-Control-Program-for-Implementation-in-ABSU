function [t,y]=RK2(tspan,y0)
%% Solving system of Differential Equations using RK2 (Second Order Runge Kutta)

a = tspan(1);
b = tspan(2);
n = 2000;     % NO. Divisions
h = (b-a)/n; % Step

y(1,:) = y0; % Initial conditions
i = 0;
for t = a:h:b
    i = i + 1;
    tt(i) = t;  % X History
    % RK2 Method
    k1 = RHS(t,y(i,:)); size(k1); size(y(i,:));
    k2 = RHS(t+h,y(i,:)+h*k1');
    phi = 0.5*k1 + 0.5*k2;
    y(i+1,:) = y(i,:) + h * phi';
end
t=tt;
