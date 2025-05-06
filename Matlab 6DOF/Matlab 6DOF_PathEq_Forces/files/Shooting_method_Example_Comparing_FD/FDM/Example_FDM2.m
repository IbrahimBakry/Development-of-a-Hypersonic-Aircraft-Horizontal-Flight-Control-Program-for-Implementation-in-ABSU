
clc; clear all

h=0.1;
t = 1:h:2.5;
N = length(t)-1;

% initial conditions
y(1) = 0.5;

n = 1; % size of the matrix
c = 1; % Number of conditions
% y(2) = 0; 
y(1) = 0.5;
for k = 1:N
if k == 1
S_k1 = (-1+h*(y(k)));
S_k  = (1-h*(y(k)));
SM   = 1;

M = k;  % Number of zeros in first and final row
% Matrix A
A = [S_k1 S_k;
     0    SM];

En = -y(k)-0.5*h*(y(k))^2;
EM = -(y(k) -2);

E(k,1) = En;
E(2*k,1) = EM; 

S = A\E;
y(k+1) = S(2); y(k) = S(1);

end
if k > 1
S_k1 = -1+h*(y(k)-y(k-1));
S_k  = 1-h*(y(k)-y(k-1));
SM = 1;

M = k;  % Number of zeros in first and final row
% Matrix A
A(k+1,k+1) = 0;
A(end,end) = SM;
A(k,k:k+1) = [S_k1 S_k];

En = -(y(k)-y(k-1)-0.5*h*(y(k)-y(k-1))^2);
EM = -(y(k) -2);

E(k+1,1)=0;
E(k,1)  = En;
E(k+1,1)= EM; 

% Solving the equations
S = A\E
y(k+1) = y(k)  + 0.5*S(k+1); y(k) = y(k-1) + 0.5*S(k); 
end
end
plot(t,y); grid on; title('Y1');  hold all