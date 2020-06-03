
close all
clear
FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
data = [];
view_chn = 11;
array_spacing = 0.75;
c_0 = 1435;
w = 64;
r = 0.01;
pthres = 0.00001;
tbe_list = [];
pe_list = [];

bandfilt = designfilt('bandpassfir','FilterOrder',1000,'CutoffFrequency1',40,'CutoffFrequency2',950,'SampleRate',12000);
set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = 100;
last_file =1500;

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);

for index = 0:(last_file-1400)/num_files-1
    
    disp([num2str(index),' / ', num2str((last_file-1400)/num_files-1)])
    
    first_file = 1400 + index*num_files;

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
    
    [peak_event_mat,loc_event_mat,num_event] = np_eSelect(data_fil,FS,w,r,pthres,view_chn,array_spacing,c_0,data_name);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/np_eSelect_plot/ACO',num2str(first_file),'_allchn.fig']]);
    %saveas(gcf,[pwd ['/np_eSelect_plot/ACO',num2str(first_file),'_chn11.png']]);
    

    loc_event = loc_event_mat(view_chn,:);
    loc_event = loc_event(~isnan(loc_event));
    peak_event = peak_event_mat(view_chn,:);
    peak_event = peak_event(~isnan(peak_event));
    time_betw_events = zeros(length(loc_event)-1,1); %time between each event

    for t = 2:length(loc_event)
        time_betw_events(t) = loc_event(t)-loc_event(t-1);
    end

    tbe_list = [tbe_list;time_betw_events];
    pe_list = [pe_list,peak_event];
    
    clear time_betw_events;
    clear peak_event;
end
