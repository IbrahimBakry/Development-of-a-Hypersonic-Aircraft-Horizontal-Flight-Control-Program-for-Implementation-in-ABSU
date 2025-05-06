function dydx = twoode(x,z)

dydx(1,1) = z(2)^2;
dydx(2,1) = z(1)^2;

% dydx(1,1) = 4*z(1)*z(2)-z(1)^2;
% dydx(2,1) = z(1)^2 + z(2)^2;

% dydx(3,1) = z(3)-z(3)*sin(z(2))*sin(z(1));
% dydx(4,1) = z(4)-z(2)*sin(z(1));

% dydx = 2*z^2;


end