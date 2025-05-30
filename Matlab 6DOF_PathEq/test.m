
syms alfav
f=@(alfav)-0.01928254697 + 0.1226892089*alfav - 0.9778696076*alfav^2;
f_dot = 0.1226916753 - 1.955703572*alfav;

N = 10000; % ? of Iterations
err = 0.001; x0 = -1;
% alfo = NRaphson(f,f_dot,N,err,x0);
[x2] = MSecant(f,N,err,x0)