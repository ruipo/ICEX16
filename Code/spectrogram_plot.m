
clear
FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
data = [];
chn = 16;
window = hamming(4096);
NFFT = length(window);

set(0,'DefaultFigureVisible','on')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
begin_file = 4565;
num_files = 60;
last_file = 4625;

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);

for index = 0:(last_file-begin_file)/num_files-1
    close all
    disp([num2str(index),' / ', num2str((last_file-begin_file)/num_files-1)])
    
    first_file = begin_file + index*num_files;

    % Start looping over ACO*.DAT files
    counter=0;
    for i = first_file:num_files+first_file

        counter=counter+1;
        filename = [prefix directory(i).name];
        fid = fopen (filename, 'r', 'ieee-le');

        if (fid <= 0)
            continue;
        end

        % Read the single precision float acoustic data samples (in uPa)
        for j = 1:NUM_CHANNELS
            aco_in(((counter-1)*NUM_SAMPLES+1):(counter*NUM_SAMPLES),j) = fread(fid, NUM_SAMPLES, 'float32');
        end

        fclose (fid);
    end

    % Nomalized to zero mean; take middle 22 channels
    data = zeros(length(aco_in),32);
    for j = 1:32
        
        data(:,j) = ((aco_in(:,j)-mean(aco_in(:,j))));
        
    end
    
    %timestamp = 1457848722.58 + first_file*2;
    timestamp = 1457914993.67 + first_file*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
    figure
    spectrogram(data(:,chn),window,[],NFFT,FS,'yaxis');
    title(['Timestamp: ',data_name,'; Channel = ',num2str(chn)]);
    set(gca,'Fontsize',20);
    caxis([50 80]);
    ylim([0.03 2.048]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/spectrogram_plot/ACO',num2str(first_file),'.fig']]);
    %saveas(gcf,[pwd ['/spectrogram_plot/ACO',num2str(first_file),'.png']]);
    
end

