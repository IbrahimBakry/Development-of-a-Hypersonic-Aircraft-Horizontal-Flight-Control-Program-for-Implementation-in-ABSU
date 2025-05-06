% initial conditions
% run before the simulink

arc     = pi/180.;

Vtot    = 1890;
chi    = 0*arc;
gamma   = 0*arc;
hn      = 30000;
phin    = 55*arc;
thetan  = 33*arc;
mf      = 81646.63;

psiV = 0.0001;   psigamma = 0.0001;   psichi = 0.0001; 
psih = 0.0001; 
psiphi = 0.0001; psitheta = 0.0001;   psim = 0.0001; 

ad = 0.1031; bd = 4.2706;
al = 0.8671; bl = 7.1112;

dynamic0 = [Vtot gamma chi hn phin thetan mf psiV psigamma...
    psichi psih psiphi psitheta psim]';


