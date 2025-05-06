function x = GSeidelm(A,B,N,err)
% Gauss-Seidel Method for Solving Linear Systems

% system we interested to solve
% A = [3    -0.1  -0.2;
%      0.1   7    -0.3;
%      0.3  -0.2   10];
% B = [7.85 ; -19.3 ; 71.4];
[m,n] = size(A);
% 
% N = 100; % Number of Iterations
% err = 0.0001; % Result Accuracy

x = zeros(14,1); % Initial value
xx(1,:) = x; % X History

for k = 2 : N
for i = 1 : n
    s = 0;
    for j = 1:n
        if j ~= i
            s = s + A(i,j) * x(j);
        end
    end
    x(i) = (1/A(i,i))*(B(i) - s);
end
xx(k,:) = x; kk = k;
Err = abs(max(xx(k,:)-xx(k-1,:))); if Err < err, break; end
end
x
disp(['The Roots are [' num2str(x) '], with accuracy: ' num2str(Err) ' ,NO Iterations: ' num2str(kk)])
