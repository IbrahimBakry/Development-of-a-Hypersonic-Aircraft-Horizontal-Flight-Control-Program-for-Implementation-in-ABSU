% Dynamic Equations
clc; clear all;

% Grid points
N = 500; a = 0; b = 100;
dt=(b-a)/(N); t = a:dt:b;
Nn = length(t);

% initial conditions
y(1,1) = 1830; y(1,2) = 0; y(1,3) = 0; y(1,4) = 30000;...
y(1,5) = 0; y(1,6) = 0; y(1,7) = 80000;

n = 7; % size of the matrix
c0 = 4; % Number of initial conditions
cf = 3; % Number of final conditions
c = c0 + cf; % % Number of total conditions

for k = 1:Nn-1
if k == 1
% Calling Matrices
[S_k1,S_k,S0,SM] = inimat(k,dt,y);
[En,E0,EM] = inidynamicf(k,dt,y);

M = n*k;  % Number of zeros
% Matrix A
A(n*k-2:n*k+c0,n*k+1:n*k+c) = 0;
A(1:c0,:) = [S0 zeros(c0,M)];
A(n*k-2:n*k+c0,n*k-c+1:n*k+n) = [S_k1 S_k];
A(n*k+c0+1:n*k+c,:) = [zeros(cf,M) SM];

% Right-hand side
E(n*k-2:n*k+c0,:)=0;
E(1:c0,1)  =  E0;
E(n*k-2:n*k+c0,1) = En;
E(n*k+c0+1:n*k+c,1) = EM; 

% Solving the equations
S = A(1:13,:)\E(1:13,1);
j2 = 7;
for i = 1:7
    j2 = j2 - 1;
    y(k+1,i) = y(k,i) + dt*S(n*k+c-j2);
end

% Checking Limits
y = LIMITS(k,y);
end

if k > 1
% Calling Matrices
[S_k1,S_k,S0,SM] = deriving_dervitevs(k,dt,y);
[En,E0,EM] = dynamicf(k,dt,y);

M = n*k;  % Number of zeros
% Matrix A
A(n*k-2:n*k+c0,n*k+1:n*k+c) = 0;
A(1:c0,:) = [S0 zeros(c0,M)];
A(n*k-2:n*k+c0,n*k-c+1:n*k+n) = [S_k1 S_k];
A(n*k+c0+1:n*k+c,:) = [zeros(cf,M) SM];

% Right-hand side
E(n*k-2:n*k+c0,:)=0;
E(1:c0,1)  =  E0;
E(n*k-2:n*k+c0,1) = En;
E(n*k+c0+1:n*k+c,1) = EM; 

% Solving the equations
S = A(1:14,:)\E(1:14,1);
j1 = 14; j2 = 7;
for i = 1:7
    j1 = j1 - 1; j2 = j2 - 1;
    y(k,i) = y(k-1,i) + dt*S(n*k+c-j1);
    y(k+1,i) = y(k,i) + dt*S(n*k+c-j2);
end

% Checking Limits
y = LIMITS(k,y);

end
end
y
subplot(211); plot(t,y(:,1)); grid on; title('V');  hold all
subplot(212); plot(t,y(:,4)); grid on; title('h');  hold all