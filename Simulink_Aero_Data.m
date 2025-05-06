%% Aerodynamic Data as Function of Alpha, Beta, Throttle

CL = xlsread('CL.xlsx');
CD = xlsread('CD.xlsx');
Cm = xlsread('Cm.xlsx');
Cy = xlsread('Cy.xlsx');
Cn = xlsread('Cn.xlsx');
Cl = zeros(length(Cn(:,1)),2);
% Isp = xlsread('Isp.xlsx');
Isp = 3146.15;
Thrust = 163.29*9.81*Isp;

alpha_CL = CL(:,1);
alpha_CD = CD(:,1);
alpha_Cm = Cm(:,1);
beta_Cy = Cy(:,1);
beta_Cn = Cn(:,1);
eta = Isp(:,1);

