
clc; clear all

N = 300;
dt=(10-0)/(N);
t = 0:dt:10;
Nn = length(t);

% initial conditions
y1(1) = 1; y2(1) = 1;

n = 2; % size of the matrix
c = 2; % Number of conditions

for k = 1:Nn-1
if k == 1
% Calling Matrices
[S_k1, S_k,S0,SM,En,E0,EM] = inimat0(k,dt,y1,y2);

M = 2*k;  % Number of zeros in first and final row
% Matrix A
A = [ S0    zeros(1,M);
      S_k1         S_k;
      zeros(1,M)    SM];
% Right-hand side
E(1,1)           = E0;
E(2*k:2*k+c/2,1) = En;
E(2*k+c,1)       = EM; 
Ed(:,k) = E;
% Solving the System
S = A\E;
y1(k+1) = S(2*k+1);
y2(k+1) = S(2*k+2);
end

if k > 1
% Calling Matrices
[S_k1,S_k,S0,SM,En,E0,EM] = nmat(k,dt,y1,y2);

M = 2*k;  % Number of zeros in first and final row
% Matrix A
A(2*k+n/2:2*k+n,2*k+n/2:2*k+n) = 0;
A(1,:) = [S0 zeros(1,M)];
A(2*k+c,:) = [zeros(1,M) SM];
A(2*k:2*k+n/2,2*k-n/2:2*k+n) = [S_k1 S_k];

% Right-hand side
E(2*k:2*k+c)=0;
E(1,1)           = E0;
E(2*k:2*k+c/2,1) = En;
E(2*k+c,1)       = EM; 
Ed(2*k+2,k) = 0; Ed(:,k) = E;

% Step-Size
% Edx2 = norm(Ed(:,k))/(length(1:k)*n), Edx1 = norm(Ed(:,k-1))/(length(1:k)*n)
% if Edx2 < Edx1
%     dt = max(h,Edx2);
% else
%     dt = h;
% end
%  dt = dt;
% Solving the equations
S = A\E;
y1(k)   = y1(k-1) + dt*S(2*k-1); 
y2(k)   = y2(k-1) + dt*S(2*k);
y1(k+1) = y1(k)   + dt*S(2*k+1); 
y2(k+1) = y2(k)   + dt*S(2*k+2);
end
end
subplot(211); plot(t,y1); grid on; title('YY1');  hold all
subplot(212); plot(t,y2); grid on; title('YY2');  hold all