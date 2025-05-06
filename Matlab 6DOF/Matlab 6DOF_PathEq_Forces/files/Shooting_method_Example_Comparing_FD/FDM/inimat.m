function [S_k1, S_k,S0,SM,En,E0,EM] = inimat(k,dt,y1,y2)
% Matrices
S_k1 = [ -1     0 ;
          0    -1];
S_k = [      1          -2*dt*(y2(k));
        -2*dt*(y1(k))         1     ];
S0 = [1 0];
SM = [0 1];

% Right-Hand Side
En = -[y1(k)-dt*(y2(k))^2;
       y2(k)-dt*(y1(k))^2];
E0 = -(y1(1) -1);
EM = -(y2(k) -2);