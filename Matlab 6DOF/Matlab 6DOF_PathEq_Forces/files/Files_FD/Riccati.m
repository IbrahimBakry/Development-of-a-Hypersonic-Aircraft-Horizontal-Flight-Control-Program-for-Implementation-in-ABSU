function X=Riccati(A,G,Q)
%RICCATI  Solves an algebraic Riccati equation
%   X = Riccati(A,G,Q) solves the algebraic Riccati equation of the form:
%       A'*X + X*A' - X*G*X + Q = 0, where X is symmetric.
n=size(A,1);
Z=[A -G
  -Q -A'];
[U1,S1]=schur(Z);
[U,S]=ordschur(U1,S1,'lhp');
X=U(n+1:end,1:n)*U(1:n,1:n)^-1;