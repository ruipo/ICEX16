
FS = 12000; 
window_size = 4096;
window = hanning(window_size);
NFFT = window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 1;
data_name = 'Hour 1';

% Set Path to DATA
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
directory = dir([prefix 'ACO0000*.DAT']);

firstfile = [600 1200 1900 1450 3050 3650 4175 4850];

for iter = 1
    iter
% first_file = 2000 + (1800*iter - 1800);
% last_file = first_file+1800;

first_file = 4850;
last_file = first_file+600;

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
% take selected channel; conver to Pa
data = aco_in(:,16)./1E6;
%data2 = [data(1:10000000,:);data(14000000:24000000,:);data(28000000:38500000,:)];
dt = 1/FS;
t = 0:dt:length(data)/FS-dt;

[psd_all,psd_med(iter,:),freq] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);
psd_75 = quantile(psd_all,0.75);
psd_25 = quantile(psd_all,0.25);
psd_90 = quantile(psd_all,0.90);
psd_10 = quantile(psd_all,0.10);

figure
hold on
plot(freq,10*log10(psd_med(iter,:)/(1E-6)^2),'linewidth',1.5)
 set(gca,'Fontsize',30);
 title([data_name,' Power Spectral Density Estimate']);
 xlabel('Frequency (Hz)');
 xlim([6 6000]);
 ylabel('Power/Frequency (dB re 1\muPa^2/Hz)');
 grid on
 hold on
plot(freq,10*log10(psd_75./1E-12),'r--','linewidth',1.5);
plot(freq,10*log10(psd_25./1E-12),'r--','linewidth',1.5);
plot(freq,10*log10(psd_90./1E-12),'k-.','linewidth',1.5);
plot(freq,10*log10(psd_10./1E-12),'k-.','linewidth',1.5); 
legend('Median','75th percentile','25th percentile','90th percentile','10th percentile')
end
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

third_oband(1,:,iter) = squeeze(mean(psd_all(:,ind1:ind2),2)); %20Hz
third_oband(2,:,iter) = squeeze(mean(psd_all(:,ind2:ind3),2)); %25Hz  
third_oband(3,:,iter) = squeeze(mean(psd_all(:,ind3:ind4),2)); %31Hz
third_oband(4,:,iter) = squeeze(mean(psd_all(:,ind4:ind5),2)); %40Hz
third_oband(5,:,iter) = squeeze(mean(psd_all(:,ind5:ind6),2)); %50Hz
third_oband(6,:,iter) = squeeze(mean(psd_all(:,ind6:ind7),2)); %63Hz
third_oband(7,:,iter) = squeeze(mean(psd_all(:,ind7:ind8),2)); %80Hz
third_oband(8,:,iter) = squeeze(mean(psd_all(:,ind8:ind9),2)); %100Hz
third_oband(9,:,iter) = squeeze(mean(psd_all(:,ind9:ind10),2)); %125Hz
third_oband(10,:,iter) = squeeze(mean(psd_all(:,ind10:ind11),2)); %160Hz
third_oband(11,:,iter) = squeeze(mean(psd_all(:,ind11:ind12),2)); %200Hz
third_oband(12,:,iter) = squeeze(mean(psd_all(:,ind12:ind13),2)); %250Hz
third_oband(13,:,iter) = squeeze(mean(psd_all(:,ind13:ind14),2)); %315Hz
third_oband(14,:,iter) = squeeze(mean(psd_all(:,ind14:ind15),2)); %400Hz
third_oband(15,:,iter) = squeeze(mean(psd_all(:,ind15:ind16),2)); %500Hz
third_oband(16,:,iter) = squeeze(mean(psd_all(:,ind16:ind17),2)); %630Hz
third_oband(17,:,iter) = squeeze(mean(psd_all(:,ind17:ind18),2)); %800Hz
third_oband(18,:,iter) = squeeze(mean(psd_all(:,ind18:ind19),2)); %1000Hz
third_oband(19,:,iter) = squeeze(mean(psd_all(:,ind19:ind20),2)); %1250Hz
third_oband(20,:,iter) = squeeze(mean(psd_all(:,ind20:ind21),2)); %1600Hz



%% Plotting

flist = [20 25 31 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600];
for j = 1:20
figure
for i = 1:8
h1 = histfit(10*log10(third_oband(j,:,i)/(1E-6)^2),[],'kernel');
x = h1(2).XData;
y = (h1(2).YData/norm(h1(2).YData,1));
plot(x,y,'linewidth',1.5);
delete(h1(1));
delete(h1(2));
hold on
end
xlabel('Power/Frequency (dB re 1\muPa^2/Hz)');
ylabel('Probability Density');
title([num2str(flist(j)),' Hz']);
set(gca,'Fontsize',30);
grid on
legend('Hour 1','Hour 2','Hour 3','Hour 4','Hour 5','Hour 6','Hour 7','Hour 8')
end

%%
for iter = 1:8
for i = 1:20
kurtlist(i) = kurtosis(10*log10(third_oband(i,:,iter)/(1E-6)^2));
end

plot(flist,kurtlist,'linewidth',2)
xlabel('Frequency (Hz)');
ylabel('Kurtosis');
set(gca,'Fontsize',30);
title(data_name);
grid on
hold on
end
%%
w = 4096;
for k = 16
    figure
     %Spectrogram usage: spectrogram(signal, window, number of overlaps,
     %number of frequency points, sampling frequency)
    spectrogram(1E6*data,hanning(w),w/2,w,12000,'yaxis');
    %set(gcf,'Position',[1          25        1600         962]);
    title(['Element ' num2str(k)])
    hcolor=colorbar;
    set(hcolor,'FontSize',30);
   
end