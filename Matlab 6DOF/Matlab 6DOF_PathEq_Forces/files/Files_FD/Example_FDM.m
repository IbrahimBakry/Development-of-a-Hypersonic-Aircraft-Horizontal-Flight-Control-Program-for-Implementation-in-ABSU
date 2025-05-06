
% 
clc; clear all

t = 1:0.5:5;
y1(1) = 1; y2(1) = 0;
n = 2; % size of the matrix
c = 2; % Number of conditions
N = length(t)-2;
for k = 1:N
if k ==1
h = t(k);
S_k1 = [   -1-h*(-0.25*cos((y1(k))/2)*(y2(k)))         -h*(-0.5*sin((y1(k))/2)+0.5*(-y2(k)))  ;
        -1-h*(-0.5*cos((y2(k))/2)*cos((y1(k))/2))  -h*(-0.5+0.5*sin((y2(k))/2)*sin((y1(k))/2))];
S_k = [   1-h*(0.25*cos((y1(k))/2)*(y2(k)))    -h*(0.25*sin((y1(k))/2)*(y2(k))+(-y2(k))/2);
       1-h*(0.5*cos((y2(k))/2)*cos((y1(k))/2))  -h*(0.5-0.5*sin((y2(k))/2)*sin((y1(k))/2))];
S0 = [1 0];
SM = [0 1]; 
  
M = 2*k;  % Number of zeros in first and final row
% Matrix A
  A = [ zeros(1,M)    S0 ;
        S_k1          S_k;
        zeros(1,M)    SM];
% Right-hand side
 En = [-h*(sin(0.5*(y1(k)))*(0.5*(y2(k)))+0.25*(y2(k))^2);
       -h*(0.5*(y2(k))+cos(0.5*(y2(k)))*sin(0.5*(y1(k))))];
 E0 = y1(1) -1;
 EM = y2(k) -3;
 
 E(1,1)     = E0;
 E(2*k:2*k+c/2,1) = En;
 E(2*k+c)   = EM; 
 
 S = A\E
 y1(k) = S(1); y2(k) = S(2); y1(k+1) = S(3); y2(k+1) = S(4); 
 
end
if k > 1
h = t(k) - t(k-1);
S_k1 = [   -1-h*(-0.25*cos((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1)))         -h*(-0.5*sin((y1(k)-y1(k-1))/2)+0.5*(y2(k-1)-y2(k)))  ;
        -1-h*(-0.5*cos((y2(k)-y2(k-1))/2)*cos((y1(k)-y1(k-1))/2))  -h*(-0.5+0.5*sin((y2(k)-y2(k-1))/2)*sin((y1(k)-y1(k-1))/2))];
S_k = [   1-h*(0.25*cos((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1)))    -h*(0.25*sin((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1))+(y2(k-1)-y2(k))/2);
       1-h*(0.5*cos((y2(k)-y2(k-1))/2)*cos((y1(k)-y1(k-1))/2))  -h*(0.5-0.5*sin((y2(k)-y2(k-1))/2)*sin((y1(k)-y1(k-1))/2))];
S0 = [1 0];
SM = [0 3];

M = 2*k;  % Number of zeros in first and final row

% Matrix A
    A(2*k+n/2:2*k+n,2*k+n/2:2*k+n) = 0;
    Finalrow = 2*k+c;
    A(1,:) = [zeros(1,M) S0];
    A(Finalrow,:) = [zeros(1,M) SM];
    A(2*k:2*k+n/2,2*k-n/2:2*k+n) = [S_k1 S_k];

% Right-hand side
 En = [-h*(sin(0.5*(y1(k)-y1(k-1)))*(0.5*(y2(k)-y2(k-1)))+0.25*(y2(k)-y2(k-1))^2);
       -h*(0.5*(y2(k)-y2(k-1))+cos(0.5*(y2(k)-y2(k-1)))*sin(0.5*(y1(k)-y1(k-1))))];
 E0 = -1;
 EM = -3;
 
 E(2*k:2*k+c)=0;
 E(1,1)     = E0;
 E(2*k:2*k+c/2,1) = En;
 E(2*k+c,1)   = EM; 
 
 S = A\E
%  y1(k) = S(2*k-1); y2(k) = S(2*k); y1(k+1) = S(2*k+1); y2(k+1) = S(2*k+2); 
 y1(k) = S(1); y2(k) = S(2); y1(k+1) = S(3); y2(k+1) = S(4);

end
end

subplot(211); plot(t(1:end-1),y1); grid on
 hold all
subplot(212); plot(t(1:end-1),y2); grid on
 hold all
