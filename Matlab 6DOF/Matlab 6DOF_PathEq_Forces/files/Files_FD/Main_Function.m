
clc; clear all

% ********************************************************
% Initial Conitions
% ********************************************************

Initial_conitions; arc = pi/180;
y0 = [1890 0.1*arc 0.1*arc 30000 55*arc 33*arc 81646.63 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001 1 1 1]';
hf = 30000; phif = 40.7*arc; thetaf = -73*arc;

%************************************************************
% Newton-Raphson
%************************************************************
n=1000; er=1; Error=0.001;
 
GO = y0(1:3);
while er >= Error
    V00 = GO(1);
    gamma00 = GO(2);
    chi00 = GO(3);
    
    [ETA,yy] = RK4(n,V00,gamma00,chi00);
    
    V = yy(n,1); gamma = yy(n,2); chi = yy(n,3);
    p1 = yy(n,15); p2 = yy(n,16); p3 = yy(n,17);
    h = yy(n,4); phi = yy(n,5); theta = yy(n,6);
    
%     J=[                      0                                         -V*sin(p2)                                         0                      ;
%               V*sin(p2)*sin(p3)/(p1+rearth)^2                   V*cos(p2)*sin(p3)/(p1+rearth)               V*sin(p2)*cos(p3)/(p1+rearth)        ;
%        -V*sin(p2)*cos(p3)/(cos(phi)*(p1+rearth)^2)      -V*cos(p2)*cos(p3)/(cos(phi)*(p1+rearth))    V*sin(p2)*sin(p3)/(cos(phi)*(p1+rearth))];
    
%       J =[1                  0                                         0                      ;
%           0     V*cos(p2)*cos(chi)/(h+rearth)                       0                         ;
%           0                  0                    -V*cos(gamma)*sin(p3)/((h+rearth)*cos(phi))];

      G=[h     -  hf;
       phi   -  phif;
       theta -  thetaf];
    GO = GO - inv(J)*G
    
    er = max(max(abs(h - hf),abs(phi - phif)),abs(theta - thetaf))
end
% disp({'f"0 = ', yy(1,3)})
%  plot(ETA,yy(:,3))
%  grid on


