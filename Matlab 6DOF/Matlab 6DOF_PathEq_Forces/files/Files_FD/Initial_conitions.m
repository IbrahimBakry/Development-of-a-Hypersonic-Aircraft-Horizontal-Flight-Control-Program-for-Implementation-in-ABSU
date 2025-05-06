% initial conditions

arc     = pi/180.;

V0    = 1890;
chi0    =  0.1*arc;
gamma0   = 0.1*arc;
h0      = 30000;
phi0    = 55*arc;
theta0  = 33*arc;
mf0      = 81646.63;

psiV0 = 0.001;   psigamma0 = 0.001;   psichi0 = 0.001; 
psih0 = 0.001; 
psiphi0 = 0.001; psitheta0 = 0.001;   psim0 = 0.001; 

p10 = 1; p20 = 1; p30 = 1;

ad = 0.1031; bd = 4.2706;
al = 0.8671; bl = 7.1112;

dynamic0 = [1890 0.1*arc 0.1*arc 30000 55*arc 33*arc 81646.63 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001 1 1 1]';

m0 = 136077.7;
m      = 136077.7;% Takeoff weight;
Sref   = 21.672;
Lref   = 8.4088;
bref   = 24.38/2;
Ac     = 27.87;
a      = 305;
rhos   = 1.225;
gs     = 9.80665;
bet    = 0.000140845;
omegae = 0.00007292;
rearth = 6.38E6;
V00    = 7606.28;
