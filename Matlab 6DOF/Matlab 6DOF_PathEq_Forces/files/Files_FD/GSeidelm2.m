function x = GSeidelm2(M,consts,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *Author: Ganindu Nanayakkara*                               %%
% *Asian Institute of Technology*                             %%
% ganindu@gmail.com / facebook.com/ganindu                    %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $$ Free For All $$                                         %%
% *Gauss Seidel Method Example*                               %%
%@paramters: Input Functions: M, x                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
% %% matrix of the problem
% M = [3 4 5;
%     4 5 3;
%     5 6 5];
%% setup
% consts = [49; 51; 61];
n = length(consts);
syms x1 x2 x3
%x_orig = [x1 x2 x3];
% x =  [1 2 5];           % enter guess here
x = x0;
%u = x;
errorc = ones(n,1);
check = zeros(n);
diag_M = diag(M);
diff_M = M - diag(diag_M);
%% check if all diagonal elements are non zero:
isDiagNonZero = all(diag_M);
if isDiagNonZero
    disp(' all diagonal elements are non zero');
else
    disp('Warning!: Non zero diagonal elements are present');
end
%% check if diagonally dominant
isDiagDominant = all(abs(diag_M - sum(diff_M,2)) >= 0);
if isDiagDominant
    disp(' Matrix Diagonally dominant');
else
    disp('Warning!: Matrix is not diagonally dominent, convergence is not guranteed');
end
%% check for symmetry
if(isequal(M,M'))
    isSymmetric = true;
    disp('Matrix is symmetric');
else
    isSymmetric = false;
    disp('Matrix is NOT symmetric');
end
%% check if positive definite
if(isSymmetric)
    if all(eig(M) >= 0)
        isPosDefinite = true;
        disp(' Matrix is Positive Definite');
    else
        isPosDefinite = false;
        disp(' Matrix is NOT positive definite');
    end
else
    disp(' Matrix is NOT positive definite');
end
%% iteration
iterations  = 0;
err = [1 1 1];
maxit = 10;
while (max(errorc) > 0.001) && (iterations < maxit)
    iterations = iterations + 1;
    ref_x = x;
    for i = 1: n
        x(i) = (consts(i) - sum(diff_M(i, :) .* x)) / M(i,i);
    end
    
    err = abs((x - ref_x) ./x ) * 100;
    err_log(:, iterations) = err;
    %      disp(['err: ', num2str(err)]);
    sol_log(:, iterations) = x;
end
log = [1:iterations; sol_log; max(err_log)]';
disp(log);
