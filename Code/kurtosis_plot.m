
close all
clear
FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
data = [];
chn = 11;
window_size = 512;
overlap = 0.5;

bandfilt = designfilt('bandpassfir','FilterOrder',1000,'CutoffFrequency1',40,'CutoffFrequency2',900,'SampleRate',12000);
set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = 15;
last_file = 2500;

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);

for index = 0:(last_file-1000)/num_files-1
    
    disp([num2str(index),' / ', num2str((last_file-1000)/num_files-1)])
    
    first_file = 1000 + index*num_files;

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

    % Nomalized to zero mean; take middle channel
    data = ((aco_in(:,chn+5)-mean(aco_in(:,chn+5))))./10^6;
    band = filtfilt(bandfilt,data);
    
    timestamp = 1457848722.58 + first_file*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
    
    kurt_ex(band,window_size,overlap,1,FS,data_name)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/psd_time_plot/ACO',num2str(first_file),'_band1.fig']]);
    saveas(gcf,[pwd ['/kurtosis_plot/ACO',num2str(first_file),'.png']]);
    
end

