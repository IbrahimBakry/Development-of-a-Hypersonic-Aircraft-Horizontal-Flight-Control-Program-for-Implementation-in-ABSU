
lat = result(:,8);
lon = result(:,7);
time = 0:255; % assume GPS measurements logged every minute

% Create a map of the route:
figure('position',[100 50 560 800])
subplot(3,1,1)
usamap([min(lat)-.05 max(lat)+.05],[min(lon)-.08 max(lon)+.08])
plotm(lat,lon,'ko-')
textm(lat(1),lon(1),'start','color',[.01 .7 .1])
textm(lat(end),lon(end),'end','color','r')