
clc; clear all

N = 30;
h=(2-1)/N;
t = 1:h:2;
Nn = length(t)-1;

% initial conditions
y1(1) = 1; y2(1) = 0;

n = 2; % size of the matrix
c = 2; % Number of conditions

T1 = 1; T2 = 2; % for changing place of Sk Sk-1

for k = 1:Nn
if k == 1
% S_k1 = [   -1-h*(-0.25*cos((y1(k))/2)*(y2(k)))         -h*(-0.5*sin((y1(k))/2)+0.5*(-y2(k)))  ;
%          -h*(-0.5*cos((y2(k))/2)*cos((y1(k))/2))  -1-h*(-0.5+0.5*sin((y2(k))/2)*sin((y1(k))/2))];
% S_k = [   1-h*(0.25*cos((y1(k))/2)*(y2(k)))    -h*(0.25*sin((y1(k))/2)*(y2(k))+(y2(k))/2);
%        -h*(0.5*cos((y2(k))/2)*cos((y1(k))/2))  1-h*(0.5-0.5*sin((y2(k))/2)*sin((y1(k))/2))];
% S0 = [1 0];
% SM = [0 1]; 
S_k1 = [-1-h*(-(y2(k))+0.5*(y1(k)))      h*(y1(k));
               h*(0.5*(y1(k)))      -1+h*(0.5*(y2(k)))];
S_k = [1-h*((y2(k))+0.5*(y1(k)))     -h*(y1(k));
            -h*(0.5*(y1(k)))       1-h*(0.5*(y2(k)))];
S0 = [1 0];
SM = [0 1];

M = 2*k;  % Number of zeros in first and final row
% Matrix A
A = [ S0 zeros(1,M);
      S_k1         S_k;
      zeros(1,M)    SM];
% Right-hand side
%  En = -[y1(k)-h*(sin(0.5*(y1(k)))*(0.5*(y2(k)))+0.25*(y2(k))^2);
%         y2(k)-h*(0.5*(y2(k))+cos(0.5*(y2(k)))*sin(0.5*(y1(k))))];
%  E0 = -(y1(1) -1);
%  EM = -(y2(k) -3);
En = -[y1(k)-h*((y1(k))*(y2(k))-0.25*(y1(k))^2);
       y2(k)-h*(0.25*(y1(k))^2 + 0.25*(y2(k))^2)];
E0 = -(y1(1) -1);
EM = -(y2(k) -3);

E(1,1)     = E0;
E(2*k:2*k+c/2,1) = En;
E(2*k+c,1) = EM; 

S = A\E;
y1(k+1) = S(3); y2(k+1) = S(4); %y1(k) = S(3); y2(k) = S(4); 

end
if k > 1
% S_k1 = [   -1-h*(-0.25*cos((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1)))         -h*(-0.5*sin((y1(k)-y1(k-1))/2)+0.5*(y2(k-1)-y2(k)))  ;
%          -h*(-0.5*cos((y2(k)-y2(k-1))/2)*cos((y1(k)-y1(k-1))/2))  -1-h*(-0.5+0.5*sin((y2(k)-y2(k-1))/2)*sin((y1(k)-y1(k-1))/2))];
% S_k = [   1-h*(0.25*cos((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1)))     -h*(0.25*sin((y1(k)-y1(k-1))/2)*(y2(k)-y2(k-1))+(-y2(k-1)+y2(k))/2);
%         -h*(0.5*cos((y2(k)-y2(k-1))/2)*cos((y1(k)-y1(k-1))/2))  1-h*(0.5-0.5*sin((y2(k)-y2(k-1))/2)*sin((y1(k)-y1(k-1))/2))];
% S0 = [1 0];
% SM = [0 1];
S_k1 = [-1-h*(-(y2(k)-y2(k-1))+0.5*(y1(k)-y1(k-1)))      h*(y1(k)-y1(k-1));
               h*(0.5*(y1(k)-y1(k-1)))              -1+h*(0.5*(y2(k)-y2(k-1)))];
S_k = [1-h*((y2(k)-y2(k-1))+0.5*(y1(k)-y1(k-1)))      -h*(y1(k)-y1(k-1));
            -h*(0.5*(y1(k)-y1(k-1)))              1-h*(0.5*(y2(k)-y2(k-1)))];
S0 = [1 0];
SM = [0 1];
M = 2*k;  % Number of zeros in first and final row
% Matrix A
A(2*k+n/2:2*k+n,2*k+n/2:2*k+n) = 0;
Finalrow = 2*k+c;
A(1,:) = [S0 zeros(1,M)];
A(Finalrow,:) = [zeros(1,M) SM];
A(2*k:2*k+n/2,2*k-n/2:2*k+n) = [S_k1 S_k];
    
% Right-hand side
%  En = -[y1(k)-y1(k-1)-h*(sin(0.5*(y1(k)-y1(k-1)))*(0.5*(y2(k)-y2(k-1)))+0.25*(y2(k)-y2(k-1))^2);
%         y2(k)-y2(k-1)-h*(0.5*(y2(k)-y2(k-1))+cos(0.5*(y2(k)-y2(k-1)))*sin(0.5*(y1(k)-y1(k-1))))];
%  E0 = -(y1(1) -1);
%  EM = -(y2(k) -3);
En = -[y1(k)-y1(k-1)-h*((y1(k)-y1(k-1))*(y2(k)-y2(k-1))-0.25*(y1(k)-y1(k-1))^2);
       y2(k)-y2(k-1)-h*(0.25*(y1(k)-y1(k-1))^2 + 0.25*(y2(k)-y2(k-1))^2)];
E0 = -(y1(1) -1);
EM = -(y2(k) -3);

E(2*k:2*k+c)=0;
E(1,1)      = E0;
E(2*k:2*k+c/2,1) = En;
E(2*k+c,1)       = EM; 

% Solving the equations
S = A\E
y1(k+1) = y1(k) + 0.001*S(2*k+1); y2(k+1) = y2(k) + 0.001*S(2*k+2); 
% y1(k)  = y1(k-1) + 0.5*S(2*k-1); y2(k) = y2(k-1) + 0.5*S(2*k);
% y1(k+1) = S(2*k+1); y2(k+1) = S(2*k+2); 
% y1(k)  =  S(2*k-1); y2(k) = S(2*k);
end
end

subplot(211); plot(t,y1); grid on; title('Y1');  hold all
subplot(212); plot(t,y2); grid on; title('Y2');  hold all

