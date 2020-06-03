path = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/ice_temp_data/MODISdata_03122016_03142016/';
directory = dir(path);
directory(1:2) = [];
num_files = length(directory);

year='2016';
mos='03';

time_vec = NaT(num_files,1);
temp_mat = zeros(406,271,num_files);
lat_mat = zeros(406,271,num_files);
long_mat = zeros(406,271,num_files);



for dd = 1:num_files
    dd
    hfile = dir([path '/' directory(dd).name '/*.hdf']);
    
    day = num2str(str2double(hfile.name(12:14))-60);
    hhmm = hfile.name(16:19);
    time_vec(dd) = datetime([year mos day hhmm],'InputFormat','yyyyMMddHHmm');
    
    hinfo = hdfinfo([hfile.folder '/' hfile.name]);
    
    for ss = 1:length(hinfo.Vgroup.Vgroup(2).SDS)
        if strcmp(hinfo.Vgroup.Vgroup(2).SDS(ss).Name,'Ice_Surface_Temperature')
            tempinfo = hinfo.Vgroup.Vgroup(2).SDS(ss);
            break
        end
    end
    
    temp = 0.01*double(hdfread(tempinfo))-273.15;
    temp = temp(3:5:end,3:5:end);
    temp_mat(:,:,dd) = temp(1:406,1:271);
    latinfo = hinfo.Vgroup.Vgroup(1).SDS(1);
    lat = hdfread(latinfo);
    lat_mat(:,:,dd) = lat(1:406,1:271);
    longinfo = hinfo.Vgroup.Vgroup(1).SDS(2);
    long = hdfread(longinfo);
    long_mat(:,:,dd) = long(1:406,1:271);
    
end
%%
[~,order] = sort(time_vec);
time_vec = time_vec(order);
temp_mat = temp_mat(:,:,order);
lat_mat = lat_mat(:,:,order);
long_mat = long_mat(:,:,order);

latlist = zeros(1,num_files);
longlist = zeros(1,num_files);
templist = zeros(1,num_files);
for dd = 1:num_files
    dd
    dist = sqrt((long_mat(1,1,dd)-73.26)^2 + (lat_mat(1,1,dd)--149.356)^2);
    for xx = 1:406
        for yy = 1:271
            if sqrt((long_mat(xx,yy,dd)-73.26)^2 + (lat_mat(xx,yy,dd)--149.356)^2) < dist
                dist = sqrt((long_mat(xx,yy,dd)-73.26)^2 + (lat_mat(xx,yy,dd)--149.356)^2);
                latlist(dd) = lat_mat(xx,yy,dd);
                longlist(dd) = long_mat(xx,yy,dd);
                templist(dd) = temp_mat(xx,yy,dd);
            end
        end
    end
end

            


%%
figure
for ii = 1:num_files
    h = pcolor(long_mat(:,:,ii),lat_mat(:,:,ii),temp_mat(:,:,ii));
    set(h,'Edgecolor','None');
    xlim([-152 -148])
    ylim([70 75])
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'fontsize',20)
    title(datestr(time_vec(ii)))
    colorbar
    hold on
    plot(-149.356, 73.26,'r*')
    pause
    clf
end

    
    