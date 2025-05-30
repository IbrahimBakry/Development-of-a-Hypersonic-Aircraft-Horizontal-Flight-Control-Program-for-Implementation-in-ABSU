%% Earth topography: changing the view point, plot a point
  % geographic coordinates of view point: Prague, Czech Republic
  latd_view_point= 55;%result(1,8);  %latitude (deg)
  lond_view_point=37;%result(1,5);  %longitude (deg)
  view_angle=2;
  
   model='ETOPO2_004arcmin'; load (model);
  [hc,hlab,name_png]=...
     rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
     'exaggeration_factor',13,'radius',6378,'units','km',...
     'cptcmap_pm','GMT_globe',...
     'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
     'graph_label',sprintf('Earth topography: (%s)',model),...
     'preview_figure_visible',0,...
     'azimuth',90+lond_view_point,'elevation',latd_view_point,...
     'view_angle',view_angle,...
     'window_height',650);

  % Plot and label the point
  % we need to know 'exaggeration_factor' and 'radius' to place the point
  %  on the displayed surface


for i = 1:1:length(result(:,5))
radius = result(i,9);
lond_view_point = result(i,5);
latd_view_point = result(i,8);

  exaggeration_factor=13;
  elev_view_point=interp2(lond_etopo2,latd_etopo2,elev_etopo2_km,...
                          lond_view_point,latd_view_point);
  r=radius/exaggeration_factor+elev_view_point;
  [xx,yy,zz]=sph2cart(lond_view_point/rad,latd_view_point/rad,r);
  hold on
  plot3(xx,yy,zz,'r.','MarkerSize',13)
%   text(xx,yy,zz,' Prague','fontSize',12,'color','red','FontWeight','bold')
  eval(sprintf('print -dpng -r0 %s;',name_png)); %save the figure again
end