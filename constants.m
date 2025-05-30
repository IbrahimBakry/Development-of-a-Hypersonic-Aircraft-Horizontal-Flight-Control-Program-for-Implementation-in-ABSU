function [m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,alX,CmX,CLX,CDX,b1X,CyX,CnX]=constants
% ********************************************************
% Function of Constant Values 
% ********************************************************

pitch     = load('Pitch.txt'); % Cm
yaw       = load('yaw.txt');   % Cn
lift      = load('lift.txt');  % CL
drag      = load('drag.txt');  % CD
sideforce = load('sideforce.txt'); % CY

alX = pitch(:,1);  % alpha
CmX = pitch(:,2);  % Cm
CLX = lift(:,2);   % CL
CDX = drag(:,2);   % CD

b1X = sideforce(:,1);  % Beta
CyX = sideforce(:,2);  % Cy
CnX = yaw(:,2);        % Cn

m      = 9300;         % Mass
Sref   = 21.672;       % Reference Area
Lref   = 8.4088;       % Reference Length
rhos   = 1.293;       
gs     = 9.80665;      % Earth Acceleration on Earth-Surface
bet    = 0.000140845;
omegae = 0.00007292;
rearth = 6.38E6;       % Earth Radius
V00    = 7606.28;

% ********************************************************
% End of Function Constants
% ********************************************************