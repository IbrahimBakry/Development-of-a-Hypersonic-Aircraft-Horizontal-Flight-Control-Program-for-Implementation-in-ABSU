function [z_dot]= BC_fun(x,z)
%% The Equation is
% 3y'' + 12y' +300y = 15*cos(50x)

z_dot(1,1) = 4*z(1)*z(2)-z(1)^2;
z_dot(2,1) = z(1)^2+z(2)^2;
z_dot(3,1) = 2*z(3);
end