function [m,Sref,Lref,rhos,gs,bet,omegae,rearth,V00,a1X,CmX,CLX,CDX,b1X,CyX,CnX] = constants
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =
% Function for Assignment of Constant Values
% = = = = = = = = = = = = = = == = == == == = == == == = == == == = == =

pitch = load('Pitch.txt');
yaw = load('yaw.txt');
lift = load('lift.txt');
drag = load('drag.txt');
sideforce = load('sideforce.txt');
a1X = pitch(:,1);
CmX = pitch(:,2);
CLX = lift(:,2);
CDX = drag(:,2);
b1X = sideforce(:,1);
CyX = sideforce(:,2);
CnX = yaw(:,2 );
m = 9300;
Sref = 21.672;
Lref = 8.4088;
rhos = 1.293;
gs = 9.80665;
bet = 0.000140845;
omegae = 0.00007292;
rearth = 6.38e6;
V00 = 7606.28;

% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **
% End of Function Constants
% * * * * * * * * * * * * ** * ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * ** **