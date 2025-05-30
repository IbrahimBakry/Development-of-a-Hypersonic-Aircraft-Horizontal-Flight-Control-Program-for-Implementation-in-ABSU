%% Simulink Constants

mass_0 = 136077.7; % [kg] Airplane weight
mass_f = 136077.7; % [kg] maximum fuel 
mass_e = 5.4431e+04; % minimum fuel

Sref = 557.42 ; % [m^2] Reference Area
Lref = 22.86 ;  % [m] Chord 
bref = 24.38/2; % [m] half-span

g0 = 9.806; % Gravity Acceleration
rho0 = 1.225;

Ix = 1572748.82;  % [kg.m^2]
Iy = 31590558.23; % [kg.m^2]
Iz = 32539630.8;  % [kg.m^2]
Ixz = 379629.026; % [kg.m^2]

arc = pi/180;

pqrt_0 = [0 0 0]';
h_theta_phi0 = [30000 55*arc 37*arc]'; h = 30000;
phi_theta_psi_0 = [0 0 0]';
Ub_Vb_Qb0 = [1830 0 0]'; ub=1830;
xme_0 = [0 0 0]';

gamma0 = -0.01;
chi0 = -0.01;

rho30000m = 0.1841;

Rearth = 6.38E6; % [m]



