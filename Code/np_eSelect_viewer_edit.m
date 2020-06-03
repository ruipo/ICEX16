 
f1 = 'event_number';
f2 = 'file_stamp';
f3 = 'start_time';
f4 = 'end_time';
f5 = 'event_dur';
f6 = 'num_arrivals';
f7 = 'peak_amp';
%f8 = 'elev_angle1';
%f9 = 'elev_angle2';

event_disp = struct(f1,[],f2,[],f3,[],f4,[],f5,[],f6,[],f7,[]);%,f8,[],f9,[]);

close all
clear
FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
data = [];
view_chn = 11;
array_spacing = 0.75;
c_0 = 1435;
w = 128;
r = 0.01;
pthres = 0.000001;
tbe_list = [];
pe_list = [];
psd_val = [];
window = hanning(4096);
window_size = length(window);
NFFT = 2*length(window);
overlap = 0.5;

bandfilt = designfilt('bandpassfir','FilterOrder',1000,'CutoffFrequency1',40,'CutoffFrequency2',950,'SampleRate',12000);
set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = 15;
last_file = 17005;

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);
count = 0;
for index = 0:(last_file-10600)/num_files-1
    close all
    disp([num2str(index),' / ', num2str((last_file-10600)/num_files-1)])
    
    first_file = 10600 + index*num_files;

    % Start looping over ACO*.DAT files
    counter=0;
    for i = first_file:num_files+first_file-1

        counter=counter+1;
        filename = [prefix directory(i).name];
        fid = fopen (filename, 'r', 'ieee-le');

        if (fid <= 0)
            continue;
        end

        % Read the single precision float acoustic data samples (in uPa)
        for j = 1:NUM_CHANNELS
            aco_in(((counter-1)*NUM_SAMPLES+1):(counter*NUM_SAMPLES),j) = fread (fid, NUM_SAMPLES, 'float32');
        end

        fclose (fid);
    end

    % Nomalized to zero mean; take middle 22 channels
    % Filter between 40-900Hz
    data = zeros(length(aco_in),22);
    data_fil = zeros(length(aco_in),22);
    for j = 1:22
        
        data(:,j) = ((aco_in(:,j+5)-mean(aco_in(:,j+5))))/1E6;
        data_fil(:,j) = filtfilt(bandfilt,data(:,j));
        
    end
    
    timestamp = 1457848722.58 + first_file*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
    
    % Event selection
    [peak_event_mat,loc_event_mat,num_event] = np_eSelect(data_fil,FS,w,r,pthres,view_chn,array_spacing,c_0,data_name);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/np_eSelect_plot/ACO',num2str(first_file),'_allchn.fig']]);
    saveas(gcf,[pwd ['/np_eSelect_plot/ACO',num2str(first_file),'_chn11.png']]);
    
    % Sort Event location
    loc_event = [];
    peak_event = [];
    for chn = 1:22
        loc_event_temp = loc_event_mat(chn,:)';
        loc_event_temp = loc_event_temp(~isnan(loc_event_temp));
        loc_event = [loc_event;loc_event_temp];
        peak_event_temp = peak_event_mat(chn,:)';
        peak_event_temp = peak_event_temp(~isnan(peak_event_temp));
        peak_event = [peak_event;peak_event_temp];
    end
    
    loc_peak_mat = zeros(length(loc_event),2);
    loc_peak_mat(:,1) = loc_event;
    loc_peak_mat(:,2) = peak_event;
    loc_peak_mat = sortrows(loc_peak_mat);
    
    loc_event = loc_peak_mat(:,1);
    peak_event = loc_peak_mat(:,2);
    
    time_betw_events = zeros(length(loc_event)-1,1);
    
    % Fill out event_disp
    
    %window = hanning(256);
    %NFFT = 2*length(window);
    %verlap = 0.5;
    %lev = -90:1:90;
    %z = -1*[7.875,7.125,6.375,5.625,4.875,4.125,3.375,2.625,1.875,1.125,0.375,-0.375,-1.125,-1.875,-2.625,-3.375,-4.125,-4.875,-5.625,-6.375,-7.125,-7.875];

    %psd_bg = psd(data,length(window),window,overlap,NFFT,FS,11,data_name);
    %xlim([40 6000]);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/psd_ts_fil/ACO',num2str(first_file),'_chn11.png']]);
    if ~isempty(loc_event);
        t_thres = 0.5;
        t_start = loc_event(1);
        ind = 1;
        for t = 2:length(loc_event)
            time_betw_events(t) = loc_event(t)-loc_event(t-1);

            if loc_event(t)-loc_event(t-1) > t_thres
                count = count+1;
                event_disp(count).f1 = count;
                event_disp(count).f2 = ['ACO',num2str(first_file)];
                event_disp(count).f3 = t_start;
                event_disp(count).f4 = loc_event(t-1);
                event_disp(count).f5 = loc_event(t-1)-t_start;
                event_disp(count).f6 = ((t-1)-ind+1)/22;
                event_disp(count).f7 = max(peak_event(ind:t-1));
                data_fil(floor(t_start*FS):ceil(loc_event(t-1)*FS),:) = NaN;
                data(floor(t_start*FS):ceil(loc_event(t-1)*FS),:) = NaN;
                %start_time = event_disp(count).f3;
                %end_time = event_disp(count).f4;
                %[event_disp(count).f8,event_disp(count).f9] = event_disp_add(data,data_fil,event_disp(count).f1,start_time,end_time,window,overlap,NFFT,FS,11,elev,z);
                t_start = loc_event(t);
                ind = t;
            end     

        end

        count = count+1; 
        event_disp(count).f1 = count;
        event_disp(count).f2 = ['ACO',num2str(first_file)];
        event_disp(count).f3 = t_start;
        event_disp(count).f4 = loc_event(end);
        event_disp(count).f5 = loc_event(end)-t_start;
        event_disp(count).f6 = (length(loc_event)-ind+1)/22;
        event_disp(count).f7 = max(peak_event(ind:end));
        data_fil(floor(t_start*FS):ceil(loc_event(end)*FS),:) = NaN;
        data(floor(t_start*FS):ceil(loc_event(end)*FS),:) = NaN;
        %start_time = event_disp(count).f3;
        %end_time = event_disp(count).f4;
        %[event_disp(count).f8,event_disp(count).f9] = event_disp_add(data,data_fil,event_disp(count).f1,start_time,end_time,window,overlap,NFFT,FS,11,elev,z);
    end
    % Fill out time between arrivals and amplitude of each arrival
    %tbe_list = [tbe_list;time_betw_events];
    %pe_list = [pe_list,peak_event];
    
    clear time_betw_events;
    clear peak_event;
    
    data_fil = data_fil(:,view_chn);
    data = data(:,view_chn);
    data_fil = data_fil(~isnan(data_fil));
    data = data(~isnan(data));
    [psd_val(index+1,:),f] = psd(data_fil,window_size,window,overlap,NFFT,FS,1,data_name);
end

%% Normplot dt = 30min

set(0,'DefaultFigureVisible','on')

[~,ind1] = min(abs(f-17.538));
[~,ind2] = min(abs(f-22.098));
[~,ind3] = min(abs(f-27.840));
[~,ind4] = min(abs(f-35.077));
[~,ind5] = min(abs(f-44.194));
[~,ind6] = min(abs(f-55.681));
[~,ind7] = min(abs(f-70.154));
[~,ind8] = min(abs(f-88.388));
[~,ind9] = min(abs(f-111.362));
[~,ind10] = min(abs(f-140.308));
[~,ind11] = min(abs(f-176.777));
[~,ind12] = min(abs(f-222.725));
[~,ind13] = min(abs(f-280.616));
[~,ind14] = min(abs(f-353.553));
[~,ind15] = min(abs(f-445.449));
[~,ind16] = min(abs(f-561.123));
[~,ind17] = min(abs(f-707.107));
[~,ind18] = min(abs(f-890.899));
[~,ind19] = min(abs(f-1122.462));
[~,ind20] = min(abs(f-1414.214));
[~,ind21] = min(abs(f-1781.797));
[~,ind22] = min(abs(f-2244.924));
[~,ind23] = min(abs(f-2828.427));
[~,ind24] = min(abs(f-3563.595));
[~,ind25] = min(abs(f-4489.848));
[~,ind26] = min(abs(f-5656.854));

third_oband = [];

third_oband(1,:) = squeeze(mean(psd_val(:,ind1:ind2),2)); %20Hz
third_oband(2,:) = squeeze(mean(psd_val(:,ind2:ind3),2)); %25Hz  
third_oband(3,:) = squeeze(mean(psd_val(:,ind3:ind4),2)); %31Hz
third_oband(4,:) = squeeze(mean(psd_val(:,ind4:ind5),2)); %40Hz
third_oband(5,:) = squeeze(mean(psd_val(:,ind5:ind6),2)); %50Hz
third_oband(6,:) = squeeze(mean(psd_val(:,ind6:ind7),2)); %63Hz
third_oband(7,:) = squeeze(mean(psd_val(:,ind7:ind8),2)); %80Hz
third_oband(8,:) = squeeze(mean(psd_val(:,ind8:ind9),2)); %100Hz
third_oband(9,:) = squeeze(mean(psd_val(:,ind9:ind10),2)); %125Hz
third_oband(10,:) = squeeze(mean(psd_val(:,ind10:ind11),2)); %160Hz
third_oband(11,:) = squeeze(mean(psd_val(:,ind11:ind12),2)); %200Hz
third_oband(12,:) = squeeze(mean(psd_val(:,ind12:ind13),2)); %250Hz
third_oband(13,:) = squeeze(mean(psd_val(:,ind13:ind14),2)); %315Hz
third_oband(14,:) = squeeze(mean(psd_val(:,ind14:ind15),2)); %400Hz
third_oband(15,:) = squeeze(mean(psd_val(:,ind15:ind16),2)); %500Hz
third_oband(16,:) = squeeze(mean(psd_val(:,ind16:ind17),2)); %630Hz
third_oband(17,:) = squeeze(mean(psd_val(:,ind17:ind18),2)); %800Hz
third_oband(18,:) = squeeze(mean(psd_val(:,ind18:ind19),2)); %1000Hz
third_oband(19,:) = squeeze(mean(psd_val(:,ind19:ind20),2)); %1250Hz
third_oband(20,:) = squeeze(mean(psd_val(:,ind20:ind21),2)); %1600Hz
third_oband(21,:) = squeeze(mean(psd_val(:,ind21:ind22),2)); %2000Hz
third_oband(22,:) = squeeze(mean(psd_val(:,ind22:ind23),2)); %2500Hz
third_oband(23,:) = squeeze(mean(psd_val(:,ind23:ind24),2)); %3150Hz
third_oband(24,:) = squeeze(mean(psd_val(:,ind24:ind25),2)); %4000Hz
third_oband(25,:) = squeeze(mean(psd_val(:,ind25:ind26),2)); %5000Hz

third_oband = 10*log10(third_oband/(1E-6)^2);
third_oband = third_oband';

data = [];
data(:,1) = third_oband(:,6);
data(:,2) = third_oband(:,7);
data(:,3) = third_oband(:,10);
data(:,4) = third_oband(:,13);
data(:,5) = third_oband(:,15);
data(:,6) = third_oband(:,17);
data(:,7) = third_oband(:,18);


figure
normplot(data)
legend({'50Hz','80Hz','160Hz','315Hz','500Hz','630Hz','800Hz'});
