%% Load Data
FS = 12000; 
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
N = 32;
 
% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
 
directory = dir([prefix 'ACO0000*.DAT']);

first_file = 9150;
last_file = 9450;
 
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

% Nomalized to zero mean;
aco_norm = zeros(length(aco_in),32);
for i = 1:N
    aco_norm(:,i) = ((aco_in(:,i)-mean(aco_in(:,i))))./10^6;
end

clear aco_in

timestamp = 1457848722.58 + first_file*2;
data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);

%%


FS = 12000; 
window_size = 2048;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
chn = 1;
data_name = 'Tape 23';

data = aco_norm(:,16);

dt = 1/FS;
t = 0:dt:length(data)/FS-dt;

[~,psd_med(16,:),freq] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);


 %%
[~,ind1] = min(abs(freq-17.538));
[~,ind2] = min(abs(freq-22.098));
[~,ind3] = min(abs(freq-27.840));
[~,ind4] = min(abs(freq-35.077));
[~,ind5] = min(abs(freq-44.194));
[~,ind6] = min(abs(freq-55.681));
[~,ind7] = min(abs(freq-70.154));
[~,ind8] = min(abs(freq-88.388));
[~,ind9] = min(abs(freq-111.362));
[~,ind10] = min(abs(freq-140.308));
[~,ind11] = min(abs(freq-176.777));
[~,ind12] = min(abs(freq-222.725));
[~,ind13] = min(abs(freq-280.616));
[~,ind14] = min(abs(freq-353.553));
[~,ind15] = min(abs(freq-445.449));
[~,ind16] = min(abs(freq-561.123));
[~,ind17] = min(abs(freq-707.107));
[~,ind18] = min(abs(freq-890.899));
[~,ind19] = min(abs(freq-1122.462));
[~,ind20] = min(abs(freq-1414.214));
[~,ind21] = min(abs(freq-1781.797));

third_oband(1,:) = squeeze(mean(psd_med(:,ind1:ind2),2)); %20Hz
third_oband(2,:) = squeeze(mean(psd_med(:,ind2:ind3),2)); %25Hz  
third_oband(3,:) = squeeze(mean(psd_med(:,ind3:ind4),2)); %31Hz
third_oband(4,:) = squeeze(mean(psd_med(:,ind4:ind5),2)); %40Hz
third_oband(5,:) = squeeze(mean(psd_med(:,ind5:ind6),2)); %50Hz
third_oband(6,:) = squeeze(mean(psd_med(:,ind6:ind7),2)); %63Hz
third_oband(7,:) = squeeze(mean(psd_med(:,ind7:ind8),2)); %80Hz
third_oband(8,:) = squeeze(mean(psd_med(:,ind8:ind9),2)); %100Hz
third_oband(9,:) = squeeze(mean(psd_med(:,ind9:ind10),2)); %125Hz
third_oband(10,:) = squeeze(mean(psd_med(:,ind10:ind11),2)); %160Hz
third_oband(11,:) = squeeze(mean(psd_med(:,ind11:ind12),2)); %200Hz
third_oband(12,:) = squeeze(mean(psd_med(:,ind12:ind13),2)); %250Hz
third_oband(13,:) = squeeze(mean(psd_med(:,ind13:ind14),2)); %315Hz
third_oband(14,:) = squeeze(mean(psd_med(:,ind14:ind15),2)); %400Hz
third_oband(15,:) = squeeze(mean(psd_med(:,ind15:ind16),2)); %500Hz
third_oband(16,:) = squeeze(mean(psd_med(:,ind16:ind17),2)); %630Hz
third_oband(17,:) = squeeze(mean(psd_med(:,ind17:ind18),2)); %800Hz
third_oband(18,:) = squeeze(mean(psd_med(:,ind18:ind19),2)); %1000Hz
third_oband(19,:) = squeeze(mean(psd_med(:,ind19:ind20),2)); %1250Hz
third_oband(20,:) = squeeze(mean(psd_med(:,ind20:ind21),2)); %1600Hz


%% Plotting

flist = [20 25 31 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600];
y1 = 63:25:238;
y2 = fliplr(y1);
y = [y1 y2];

for j = 1:20
figure
plot(10*log10(third_oband(j,:)/(1E-6)^2),y,'linewidth',1.5);
xlabel('Power/Frequency (dB re 1\muPa^2/Hz)');
ylabel('Depth (m)');
ylim([77 252]);
xlim([40 90]);
set (gca,'Ydir','reverse');
set(gca,'Fontsize',30);
title([num2str(flist(j)),' Hz']);
grid on

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf,[pwd ['/SPL_depth/SPL_depth_',num2str(flist(j)),'Hz.fig']]);
end
