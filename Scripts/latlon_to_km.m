usgs_distance = [];
lat_usgs = [];
lon_usgs = [];
for i = 1:length(data.ifis_id)
    R = 6371; %radius of earth
    lat_diff = (usgs_lat - data.lat(i))*pi()/180;
    lon_diff = (usgs_lon - data.lon(i))*pi()/180;
    a = sin(lat_diff/2).*sin(lat_diff/2) + cos(data.lat(i)*pi()/180).*cos(usgs_lat*pi()/180).*sin(lon_diff/2).*sin(lon_diff/2);
    c = 2*atan2(a.^(1/2),(1-a).^(1/2));
    d = R*c;
    [usgs_distance(i,1),ind] = min(d);
    lat_usgs(i,1) = usgs_lat(ind);
    lon_usgs(i,1) = usgs_lon(ind);
end