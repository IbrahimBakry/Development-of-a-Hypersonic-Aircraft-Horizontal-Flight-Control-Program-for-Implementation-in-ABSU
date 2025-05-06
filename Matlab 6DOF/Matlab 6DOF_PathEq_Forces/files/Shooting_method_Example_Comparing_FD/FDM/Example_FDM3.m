J = 10; dx = (1-0)/J; x = (0:dx:1)';
b = zeros(J+1,1);
b(2:J) = 12 * dx^2 * x(2:J).^2;
A = sparse(J+1,J+1);
A(1,1) = 1.0; A(J+1,J+1) = 1.0;
for j=2:J
A(j,[j-1, j, j+1]) = [1, -2, 1];
end
% solving the matrix problem looks like:
Y = A \ b; % solve A Y = b
% also get exact soln on fine grid:
xf = 0:1/1000:1; yexact = xf.^4 - xf;
plot(x,Y,'o','markersize',12);hold on; plot(xf,yexact)
grid on, xlabel x, legend('finite diff','exact')