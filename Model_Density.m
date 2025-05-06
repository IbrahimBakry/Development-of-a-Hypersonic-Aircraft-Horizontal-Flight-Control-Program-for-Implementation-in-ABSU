% Density Model
rhos = 1.225;
h = 0:1000:40000;
H = 10351.8-(0.0368512).*h-(1.02368e-5).*h.^2+(2.63363e-10).*h.^3;
rho1 = rhos.*exp(-h./H);

% IS0 2533
R = 287.058; i = 0;
for j=0:1000:40000
    i = i+1;
if j < 11000
    T = 288.15 - 0.0065.*j;
    P = 101325.*(T./288.15).^5.2559;
elseif j>=11000
    T = 216;
    P = 22630.*exp(-0.00015769.*(j-11000));
end
rho2(i) = P./(R.*T);
end

plot(h,rho1,'-g',h,rho2,'--ro',...
        'LineWidth',2,...
        'MarkerEdgeColor','b')
grid on
xlabel('Hight')
ylabel('Density')
axis([0 40000 -0.2 1.3 ])
legend('Full Function','ISO-2533')