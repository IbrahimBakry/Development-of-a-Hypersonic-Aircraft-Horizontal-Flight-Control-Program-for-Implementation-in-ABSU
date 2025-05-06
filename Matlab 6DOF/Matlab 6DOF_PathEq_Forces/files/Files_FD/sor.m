function [u,it] = sor(A,F,w,u0,dtol,itmax)
% SOR solves the linear system Au=F using successive over relaxation. Where
% A is diagonally dominant square matrix.
%
% syntax
%
% u = sor(A,F) finds the solution of the linear system applying successive
% over relaxation technique. In this case the parameter w, the initial guess
% u0, the stopping criterion dtol and the maximum number of iterations
% itmax are to default values: w = 1.5, u0 = F.*0+1, dtol = 1e-3, itmax =
% 1000.
%
% The parameter w can theoretically take values between 0 - 2. However, the
% number of iterations may be reduced if w is taken near 1.5. Using w = 1
% will be equivalent to use Gauss-Seidel method. 
%
% u0 is the initial guess used by SOR to find the solution. The better the
% initial guess the faster the solution is found. However, if the matrix A
% is strictly diagonally dominant, the method converges to the solution no
% matter the initial guess (see
% http://ocw.mit.edu/OcwWeb/Aeronautics-and-Astronautics/16-920JNumerical-M
% ethods-for-Partial-Differential-EquationsSpring2003/CourseHome/index.htm
% Lecture 6 "Iterative Techniques").
%
% As stated before the matrix A should be diagonally dominant. This program
% DOESN'T check this property; you will just get a NaN vector or something 
% similar as an answer.
% 
% Finally the program will stop when abs(u-ui)/abs(ui) < dtol, where u and
% ui are the solution and initial guess used for a particular iteration. 
%
% Example 1
% e = ones(5,1);
% A = spdiags([-e 2*e -e],-1:1,5,5);
% F = [1:5]';
% [u, it] = sor(A,F,[],[],1e-5)
%
% u =
%     5.8333
%    10.6666
%    13.5000
%    13.3333
%     9.1667
%     
% it =
%     18
% 
% Comparing this with gaussian elimination we get:
% 
% A\F
%
% ans =
%     5.8333
%    10.6667
%    13.5000
%    13.3333
%     9.1667
% 
% Example 2 (what happens when A is not diagonal dominant?) 
%
% [u, it] = sor(magic(5),F,[],[],1e-5)
% u =
%    NaN
%    NaN
%    NaN
%    NaN
%    NaN
% 
% it =
%    277
% ===============================================
%% Checking proper values of the input parameters
% ===============================================
if nargin < 2
    error ('Too few arguments')
end
if nargin < 6 || isempty(itmax) == 1
    itmax = 1000;
end
if nargin < 5 || isempty(dtol) == 1
    dtol = 1e-3;
end
if nargin < 4 || isempty(u0) == 1
    u0 = F*0 + 1;
end
if nargin < 3 || isempty(w) == 1
    w = 1.5;
elseif 0 > w || w > 2
    warning('improper value for the parameter w. w = 1.5 was taken instead')
    w = 1.5;
end
if size(A,1) ~= size(A,2)
    error('Matrix A should be square.')
elseif size(A,1) ~= size(F,1)
    error('Mismatch between dimensions of A and F')
end
%===================
% SOR Begins
%===================
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
DL = D - w*L;
UE = (1-w)*D + w*U;
FE = w*F;
for it = 1:itmax
    rhv = UE*u0 + FE;
    u = DL\rhv;
    if abs(u-u0)/abs(u0) < dtol
        return
    else
        u0 = u;
    end
end
warning('Maximun number of iterations reached. Results may not be accurate')