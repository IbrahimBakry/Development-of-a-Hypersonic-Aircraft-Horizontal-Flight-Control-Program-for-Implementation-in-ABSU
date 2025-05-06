function [t,y]=RK2(tspan,y0,ii)
%% Solving system of Differential Equations using RK2 (Second Order Runge Kutta)

a = tspan(1);
b = tspan(2);
n = 10e9;     % NO. Divisions
h = (b-a)/n; % Step

y(1,:) = y0; % Initial conditions
i = 0; ii = 0;
for t = a:h:b
    i = i + 1;
    tt(i) = t;  % X History
    % RK2 Method
    ii = ii +1;
    k1 = RHS(t,y(i,:),ii); size(k1); size(y(i,:));
    k2 = RHS(t+h,y(i,:)+h*k1',ii);
    phi = 0.5*k1 + 0.5*k2;
    y(i+1,:) = y(i,:) + h * phi';
end
t=tt;
t
