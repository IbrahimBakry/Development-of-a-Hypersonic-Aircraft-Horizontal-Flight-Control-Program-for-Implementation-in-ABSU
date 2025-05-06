function [ETA,zz]=RK4(n,V00,gamma00,chi00)
% Eng. Ibrahim Bakry
% For Public Use
% Done at 5/1/2016

a=0;
b=100;   
h=(b-a)/n;    % eta Pitch;
arc = pi/180;
z=[V00 gamma00 chi00 30000 55*arc 33*arc 81646.63 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001 1 1 1]';

ii=1000;
i=0;
for eta=a:h:b
    i=i+1;
    ETA(i)=eta;
    zz(i,:)=z';
    
    k1=h*dynamic2(eta,z);
    k2=h*dynamic2(eta+h/2,z+k1/2);
    k3=h*dynamic2(eta+h/2,z+k2/2);
    k4=h*dynamic2(eta+h,z+k3);
    k=1/6*(k1+2*k2+2*k3+k4);
    z=z+k;
    
    if i==ii
       i 
        ii=ii+1000;
    end
end
end

