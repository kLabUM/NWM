clear all, close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM');
lat = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','lat');
lon = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','lon');
ID = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','link');
to = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','to');

states = shaperead('usastatelo' , 'UseGeoCoords' , true);
ind = inpolygon(lon,lat,states(15).Lon,states(15).Lat);
lat = lat(ind);
lon = lon(ind);
ID = ID(ind);
to = to(ind);

% usgs = readtable('Scripts/USGS_sites.txt');
usgs = readtable('Scripts/inventory');
usgs = usgs(2:end,:);
usgs_lat = cellfun(@str2double,usgs.dec_lat_va);
usgs_lon = cellfun(@str2double,usgs.dec_long_va);
ind = inpolygon(usgs_lon,usgs_lat,states(15).Lon,states(15).Lat);
usgs_lon = usgs_lon(ind);
usgs_lat = usgs_lat(ind);

ifis = readtable('Scripts/IFIS_NWM_data_new.txt');


hold on
plot(states(15).Lon,states(15).Lat,'k')
for i = 1:length(ID)
    to_ind = find(to==ID(i));
        if isempty(to_ind)==0
            for j=1:length(to_ind)
                plot([lon(i);lon(to_ind(j))],[lat(i);lat(to_ind(j))],'b')
            end
        end
end

h1 = plot(ifis.lon,ifis.lat,'o','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k');
h2 = plot(usgs_lon,usgs_lat,'o','MarkerSize',10,'MarkerFaceColor','c','MarkerEdgeColor','k');
ind = find(strcmp(ifis.ifis_id,'MQKTARV03'));
h3 = plot(ifis.lon(ind),ifis.lat(ind),'d','MarkerSize',20,'MarkerFaceColor','g','MarkerEdgeColor','k');
ind = find(strcmp(ifis.ifis_id,'RAPIDTRB01'));
h4 = plot(ifis.lon(ind),ifis.lat(ind),'d','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor','k');
ind = find(strcmp(ifis.ifis_id,'DSMNSRV02'));
h5 = plot(ifis.lon(ind),ifis.lat(ind),'d','MarkerSize',20,'MarkerFaceColor','r','MarkerEdgeColor','k');
legend([h1,h2,h3,h4,h5],{'IFIS Gages','USGS Gages','Example Site 1','Example Site 2','Example Site 3'})
xlabel('Longitude')
ylabel('Latitude')
hold off
set(gca,'FontSize',20)
%%

file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/Figures/ifis_usgs_locations.png';
export_fig(file,'-nocrop','-transparent','-r100',gcf)