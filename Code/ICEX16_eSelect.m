
%% Load Data
FS = 12000; 
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
 
% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
 
directory = dir([prefix 'ACO0000*.DAT']);

first_file = 2000;
last_file = first_file + 180;
 
% Read DATA
aco_in = zeros(NUM_SAMPLES * (last_file-first_file), 32);
  
% Start looping over ACO*.DAT files
counter=0;
for i = first_file:last_file-1
 
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
aco_norm = zeros(length(aco_in),32);
for i = 1:NUM_CHANNELS
    aco_norm(:,i) = ((aco_in(:,i)-mean(aco_in(:,i))))./10^6;
end

clear aco_in

timestamp = 1457848722.58 + first_file*2;
data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
%%

aco = [aco_norm(1:10000000,:);aco_norm(14000000:24000000,:);aco_norm(28000000:38500000,:)];
%%
bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',80,'CutoffFrequency2',160,'SampleRate',12000);

data = filtfilt(bandfilt,aco);
data_fil = data;


FS = 12000;         
NUM_CHANNELS = 32;
view_chn = 16;
array_spacing = .75;
c_0 = 1440;
w = 10;
r = 0.01;
pthres = 0.00001;
tbe_list = [];
data_name = 'Hour 1';

f1 = 'event_number';
f2 = 'start_time';
f3 = 'end_time';
f4 = 'event_dur';
f5 = 'num_arrivals';
f6 = 'peak_amp';
event_disp = struct(f1,[],f2,[],f3,[],f4,[],f5,[],f6,[]);

[peak_event_mat,loc_event_mat,num_event] = np_eSelect(data,FS,w,r,pthres,view_chn,array_spacing,c_0,data_name);

%%
% Sort Event location
loc_event = [];
peak_event = [];
for chn = 1:32
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

count = 0;
if ~isempty(loc_event);
    t_thres = 0.05;
    t_start = loc_event(1);
    ind = 1;
    for t = 2:length(loc_event)
        time_betw_events(t-1) = loc_event(t)-loc_event(t-1);

        if loc_event(t)-loc_event(t-1) > t_thres
            count = count+1;
            event_disp(count).event_number = count;
            event_disp(count).start_time = t_start;
            event_disp(count).end_time = loc_event(t-1);
            event_disp(count).event_dur = loc_event(t-1)-t_start;
            event_disp(count).num_arrivals = ((t-1)-ind+1)/32;
            event_disp(count).peak_amp = max(peak_event(ind:t-1));
            %data(floor(t_start*FS):ceil(loc_event(t-1)*FS),:) = NaN;
            data_fil(floor(t_start*FS)+1:ceil(loc_event(t-1)*FS),:) = NaN;
            t_start = loc_event(t);
            ind = t;
        end     

    end

    count = count+1; 
    event_disp(count).event_number = count;
    event_disp(count).start_time = t_start;
    event_disp(count).end_time = loc_event(end);
    event_disp(count).event_dur = loc_event(end)-t_start;
    event_disp(count).num_arrivals = (length(loc_event)-ind+1)/32;
    event_disp(count).peak_amp = max(peak_event(ind:end));
    %data_fil(floor(t_start*FS):ceil(loc_event(end)*FS),:) = NaN;
    data_fil(floor(t_start*FS)+1:ceil(loc_event(end)*FS),:) = NaN;
end 