function [alf]=alpha_equation(V,psiV,psigamma)
% Solving Alpha Equation  

ad = 0.1031; bd = 4.2706;
al = 0.8671; bl = 7.1112;
alf = double((psiV*ad - (psigamma/V)*al)/((psigamma/V)*bl - psiV*bd))

% syms alfav
% dCL1 = 0.8671 + 7.1112*alfav;
% dCD1 = 0.1031 + 4.2706*alfav;
% 
% dCL2 =  7.1112;
% dCD2 =  4.2706;
% 
% % Alpha Equation f(alpha)
% f = - psiV*(dCD1) + psigamma*(dCL1/V);
% 
% % Equation f_dot(alpha)
% f_dot = - psiV*(dCD2) + psigamma*(dCL2/V);
% 
% N = 1000; % ? of Iterations
% err = 0.001; x0 = 0.1;
% alfo = NRaphson(f,f_dot,N,err,x0);

% arc = pi/180;
% if alfo <= (-03)*arc, alfo = -03*arc; end
% if alfo >= (+21)*arc, alfo = +21*arc; end