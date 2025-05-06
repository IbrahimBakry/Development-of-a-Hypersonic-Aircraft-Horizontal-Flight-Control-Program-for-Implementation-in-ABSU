function [xx,zz]=RK4shoot(n,alpha)
% Here we Just Solve The Two Equation System and
%  we need Just "alpha" and "n"

a=0;
b=5;
h=(b-a)/n;
z=[1;alpha;1];

i=0;
for x=a:h:b
    i=i+1;
    xx(i)=x;
    zz(i,:)=z';
    
    k1=h*BC_fun(x,z);
    k2=h*BC_fun(x+h/2,z+k1/2);
    k3=h*BC_fun(x+h/2,z+k2/2);
    k4=h*BC_fun(x+h,z+k3);
    k=(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    z = z + k;
end


