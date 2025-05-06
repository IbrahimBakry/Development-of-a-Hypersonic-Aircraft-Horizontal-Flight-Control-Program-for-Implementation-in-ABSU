function [x]=NRaphsonb(f,f_dot,N,err,x0)
%%% Newton Raphson Method for Roots of Equations
% This method is faster than all the other methods
% for solving Equations. And the most accurate.  
% N: ? Iterations
% err: % Result Accuracy
% f: Eq. want to solve 
% x0 = 0; % Initial Value;

x = x0;
xx(1) = x0; % X History 
for i = 2:N
    x = x - subs(f,x)/subs(f_dot,x); % Newton Raphson Law.
    xx(i) = x; ii = i-1;
    Err = abs(xx(i) - xx(i-1)); if Err < err, break; end % Condition of convergenc
end
disp(['beta is: ' num2str(double(x)) ' ,with accuracy: ' num2str(double(Err)) ' ,NO Iterations: ' num2str(ii)])
