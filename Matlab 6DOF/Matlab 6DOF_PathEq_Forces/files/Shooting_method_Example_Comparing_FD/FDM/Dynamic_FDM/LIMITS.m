function y = LIMITS(k,y)
% Velocity V:
if y(k+1,1)>=2440, y(k+1,1)=2440; end
if y(k+1,1)<=1525, y(k+1,1)=1525; end
% Gamma:
if y(k+1,2)>=89/57.3, y(k+1,2)=89/57.3; end
if y(k+1,2)<=-89/57.3, y(k+1,2)=-89/57.3; end
if y(k+1,2)==0, y(k+1,2)=0.01; end
% Chi:
if y(k+1,3)>=3.14, y(k+1,3)=3.14; end
if y(k+1,3)<=-3.14, y(k+1,3)=-3.14; end
if y(k+1,3)==0, y(k+1,3)=0.01; end
% Hight h:
% if y(k+1,4)>=30500, y(k+1,4)=30500; end
% if y(k+1,4)<=29500, y(k+1,4)=29500; end
if y(k+1,4)>30000 || y(k+1,4)<30000, y(k+1,4)=30000; end
% Phi:
if y(k+1,5)>=89/57.3, y(k+1,5)=89/57.3; end
if y(k+1,5)<=-89/57.3, y(k+1,5)=-89/57.3; end
if y(k+1,5)==0, y(k+1,5)=0.01; end
% Theta:
if y(k+1,6)>=3.14, y(k+1,6)=3.14; end
if y(k+1,6)<=-3.14, y(k+1,6)=-3.14; end
if y(k+1,6)==0, y(k+1,6)=0.01; end
% MF:
if y(k+1,7)>=80000, y(k+1,7)=80000; end
if y(k+1,7)<=0, y(k+1,7)=1; end


% Checking dividing on Zero:
% parameters in Denomerator are V, gamma, phi
for i = 1:7
    if y(k+1) == y(k),
        y(k+1) = y(k+1)+0.001;
    end
end




