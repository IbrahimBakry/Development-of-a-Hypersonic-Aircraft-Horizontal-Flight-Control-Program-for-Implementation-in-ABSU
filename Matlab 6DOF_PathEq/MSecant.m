function [x2]=MSecant(f,N,err,x0)
%% Modified Secant Method for Roots of Equations
%V Number of Iterations
%err Accuracy of Result;
% f  Wanted Eq. to be Solved
% x0 Initial Condition;
delta = 0.01; % small perturbation factor
xx(1) = x0; % History of x;

for i = 2 : N
    x2 = x0 - (f(x0)*delta*x0)/(f(x0+delta*x0) - f(x0)); % Modified Secant Law
    x0 = x2;
    xx(i) = x2; ii = i; % Number of Iterations
    Err = abs(xx(i) - xx(i-1)); if Err < err, break; end % Convergenc Condition
end
disp(['The Root is: ' num2str(x2) ' ,with accuracy: ' num2str(Err) ' ,No Iterations: ' num2str(ii)])
