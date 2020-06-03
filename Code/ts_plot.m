
%close all
clear
FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 11;

bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',200,'CutoffFrequency2',1000,'SampleRate',12000);

set(0,'DefaultFigureVisible','on')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/engin_test/DURIP/DURIP_20290817T014032/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
first_file = 4565;
num_files = 60;
last_file = 4625;

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);

%for index = 0:(last_file)/num_files-1
    
    %disp([num2str(index),' / ', num2str((last_file)/num_files-1)])
    
    %first_file = 1 + index*num_files;

    % Start looping over ACO*.DAT filess
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
            aco_in(((counter-1)*NUM_SAMPLES+1):(counter*NUM_SAMPLES),j) = fread (fid, NUM_SAMPLES, 'float32');
        end

        fclose (fid);
    end

    % Nomalized to zero mean; take selected channel
        
    data = ((aco_in(:,chn+5)-mean(aco_in(:,chn+5))))./10^6;
       
    timestamp = 1457848722.58 + first_file*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);

    time = (1/(FS))*(0:length(data)-1/FS);
    band1 = filtfilt(bandfilt,data);
    
    figure
    plot(time,band1);
    set(gca,'Fontsize',20);
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    title(['Filtered TS (1300-1400Hz); Timestamp = ',data_name]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/ts_plot/ACO',num2str(first_file),'00.fig']]);
    %saveas(gcf,[pwd ['/ts_plot/ACO',num2str(first_file),'.png']]);
    
%end
    
% band1filt = designfilt('bandpassfir','FilterOrder',1000,'CutoffFrequency1',40,'CutoffFrequency2',900,'SampleRate',12000);
% band2filt = designfilt('bandpassfir','FilterOrder',2000,'CutoffFrequency1',40,'CutoffFrequency2',80,'SampleRate',12000);
% band3filt = designfilt('bandpassfir','FilterOrder',2000,'CutoffFrequency1',80,'CutoffFrequency2',160,'SampleRate',12000);
% band4filt = designfilt('bandpassfir','FilterOrder',2000,'CutoffFrequency1',160,'CutoffFrequency2',320,'SampleRate',12000);
% band5filt = designfilt('bandpassfir','FilterOrder',2000,'CutoffFrequency1',320,'CutoffFrequency2',640,'SampleRate',12000);
% band6filt = designfilt('bandpassfir','FilterOrder',2000,'CutoffFrequency1',640,'CutoffFrequency2',1280,'SampleRate',12000);

% band1 = filtfilt(band1filt,data);
% band2 = filtfilt(band2filt,data);
% band3 = filtfilt(band3filt,data);
% band4 = filtfilt(band4filt,data);
% band5 = filtfilt(band5filt,data);
% band6 = filtfilt(band6filt,data);

% figure
% plot(time,band1);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered TS (40-900Hz);');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(gcf,[pwd ['/ts_plot/ACO',num2str(start),'00.fig']]);
% saveas(gcf,[pwd ['/ts_plot/ACO',num2str(start),'00.png']]);

% figure
% plot(time,band2);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered Time Series (40-80Hz)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure(2),[pwd ['/ts_plot/ACO',num2str(start),'00_band2.fig']]);
% saveas(figure(2),[pwd ['/ts_plot/ACO',num2str(start),'00_band2.png']]);
% 
% figure
% plot(time,band3);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered Time Series (80-160Hz)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure(3),[pwd ['/ts_plot/ACO',num2str(start),'00_band3.fig']]);
% saveas(figure(3),[pwd ['/ts_plot/ACO',num2str(start),'00_band3.png']]);
% 
% figure
% plot(time,band4);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered Time Series (160-320Hz)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure(4),[pwd ['/ts_plot/ACO',num2str(start),'00_band4.fig']]);
% saveas(figure(4),[pwd ['/ts_plot/ACO',num2str(start),'00_band4.png']]);
% 
% figure
% plot(time,band5);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered Time Series (320-640Hz)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure(5),[pwd ['/ts_plot/ACO',num2str(start),'00_band5.fig']]);
% saveas(figure(5),[pwd ['/ts_plot/ACO',num2str(start),'00_band5.png']]);
% 
% figure
% plot(time,band6);
% set(gca,'Fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% title('Filtered Time Series (640-1280Hz)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure(6),[pwd ['/ts_plot/ACO',num2str(start),'00_band6.fig']]);
% saveas(figure(6),[pwd ['/ts_plot/ACO',num2str(start),'00_band6.png']]);

