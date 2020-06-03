path = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/ice_temp_data/VIIRSdata_03122016_03142016/';
directory = dir(path);
directory(1:2) = [];
num_files = length(directory);

year='2016';
mos='03';

time_vec = NaT(num_files,1);
temp_mat = zeros(3200,3232,num_files);
lat_mat = zeros(3200,3232,num_files);
long_mat = zeros(3200,3232,num_files);



for dd = 1:num_files
    dd
    ncfile = dir([path '/' directory(dd).name '/*.nc']);
    
    day = num2str(str2double(ncfile.name(12:14))-60);
    hhmm = ncfile.name(16:19);
    time_vec(dd) = datetime([year mos day hhmm],'InputFormat','yyyyMMddHHmm');
    
    temp = ncread([ncfile.folder '/' ncfile.name],'/IST_Data/IST_map');
    lat = ncread([ncfile.folder '/' ncfile.name],'/Geolocation_Data/latitude');
    long = ncread([ncfile.folder '/' ncfile.name],'/Geolocation_Data/longitude');
    temp(long>0) = NaN;
    lat(long>0) = NaN;
    long(long>0) = NaN;
    
    %temp = 0.01*double(hdfread(tempinfo))-273.15;
    temp = temp(1:3200,1:3232);
    temp(temp>655) = NaN;
    temp(temp<=50) = NaN;
    temp_mat(:,:,dd) = temp;
    lat_mat(:,:,dd) = lat(1:3200,1:3232);
    long_mat(:,:,dd) = long(1:3200,1:3232);
    
end

%%
[~,order] = sort(time_vec);
time_vec = time_vec(order);
temp_mat = temp_mat(:,:,order);
lat_mat = lat_mat(:,:,order);
long_mat = long_mat(:,:,order);

figure
for ii = 11
    h = pcolor(long_mat(:,:,ii),lat_mat(:,:,ii),temp_mat(:,:,ii)-273.15);
    set(h,'Edgecolor','None');
    xlim([-150 -148])
    ylim([72 74])
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'fontsize',20)
    title(datestr(time_vec(ii)))
    caxis([-30,-15])
    colorbar
    a = colorbar;
    a.Label.String = '^o C';
    colormap jet
    hold on
    plot(-149.356, 73.26,'kp','MarkerFaceColor','black','MarkerSize',10)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,['/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/ice_temp_data/figures/' num2str(ii) '.png']);
    %pause
    %clf
end
%%
latlist = zeros(1,num_files);
longlist = zeros(1,num_files);
templist = zeros(1,num_files);
for dd = 1:num_files
    dd
    dist = sqrt((long_mat(1,1,dd)-73.26)^2 + (lat_mat(1,1,dd)--149.356)^2);
    for xx = 1:3200
        for yy = 1:3232
            if (sqrt((long_mat(xx,yy,dd)-73.26)^2 + (lat_mat(xx,yy,dd)--149.356)^2) < dist) && (temp_mat(xx,yy,dd)<=310) && (210<=temp_mat(xx,yy,dd))
                dist = sqrt((long_mat(xx,yy,dd)-73.26)^2 + (lat_mat(xx,yy,dd)--149.356)^2);
                latlist(dd) = lat_mat(xx,yy,dd);
                longlist(dd) = long_mat(xx,yy,dd);
                templist(dd) = temp_mat(xx,yy,dd);
            end
        end
    end
end