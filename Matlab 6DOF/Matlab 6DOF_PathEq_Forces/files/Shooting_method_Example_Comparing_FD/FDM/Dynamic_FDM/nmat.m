function [S_k1,S_k,S0,SM,En,E0,EM] = nmat(k,dt,y1,y2)
S_k1 = [           -1           0.5*dt*(y2(k)-y2(k-1));
        0.5*dt*(y1(k)-y1(k-1))           -1          ];
S_k = [             1          -0.5*dt*(y2(k)-y2(k-1));
        0.5*dt*(y1(k)-y1(k-1))            1          ];
S0 = [1 0];
SM = [0 1];

En = -[y1(k)-y1(k-1)-0.25*dt*(y2(k)-y2(k-1))^2;
       y2(k)-y2(k-1)-0.25*dt*(y1(k)-y1(k-1))^2];
E0 = -(y1(1) -1);
EM = -(y2(k) -2);