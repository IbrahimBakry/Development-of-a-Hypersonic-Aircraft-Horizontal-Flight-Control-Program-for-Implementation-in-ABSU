load coastlines
axesm mollweid
framem('FEdgeColor','blue','FLineWidth',0.5)

plotm(coastlat,coastlon,'LineWidth',1,'Color','blue')

% figure
% geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
% axesm balthsrt

figure;
axesm sinusoid
framem on; gridm on; tightmap tight
load coast
patchm(lat, long,'g')
citylats = result(:,8); citylongs = result(:,7);
plotm(citylats,citylongs,'r*')

