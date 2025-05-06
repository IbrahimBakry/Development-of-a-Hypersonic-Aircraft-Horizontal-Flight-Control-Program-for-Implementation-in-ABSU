clc; clear all

arc = pi/180;
x = [1890 0.1*arc 0.1*arc 30000 55*arc 33*arc 81646.63 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001]';
xd = [1890 1*arc 1*arc 30000 40.7*arc -73*arc 1000 0 0 0 0 0 0 0]';

for i = 1:10
ii(i) = i;

% x(:,i+1) = dynamic3(i,x);
V1 = x(1,i);         gamma1 = x(2,i);  chi1  = x(3,i);
h1 = x(4,i);         phi1   = x(5,i);  theta1= x(6,i);
mf = x(7,i);         psiV   = x(8,i);  psigamma= x(9,i);
psichi = x(10,i);    psih   = x(11,i); psiphi= x(12,i);
psitheta = x(13,i);  psim   = x(14,i);

[A,B] = LinAB(V1,gamma1,chi1,h1,phi1,theta1,mf,psiV,psigamma,psichi,psih...
    ,psiphi,psitheta,psim);

Q = 10*eye(14); R = 1*ones(1,1);
% [K] = lqr(A,B,Q,R);
G = B*inv(R)*B'
p=Riccati(A,G,Q)
k = inv(R)*B'*p
u = -k*x(:,i)

if i >=2, 
    u = k*(x(:,i)-xd)
end

x(:,i+1) = A*x(:,i) + B*u

end