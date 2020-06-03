
f0 = 'file_name';
f1 = 'band1_slope';
f2 = 'band1_yint';
f3 = 'band1_r2';
f4 = 'band2_slope';
f5 = 'band2_yint';
f6 = 'band2_r2';

psd_fit = struct(f0,[],f1,[],f2,[],f3,[],f4,[],f5,[],f6,[]);

FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
data = [];
chn = 11;
window = hamming(4096);
window_size = length(window);
NFFT = 2*length(window);
overlap = 0.5;
freq = 0:FS/NFFT:FS/2;

set(0,'DefaultFigureVisible','on')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = 15;
last_file = 17005;
begin_file = 1000;
psd_val = [];
tlist = [];

% Read DATA
aco_in = zeros (NUM_SAMPLES * num_files, 32);

for index = 0:(last_file-begin_file)/num_files-1
    %close all
    disp([num2str(index),' / ', num2str((last_file-begin_file)/num_files-1)])
    
    first_file = begin_file + index*num_files;

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
    data = zeros(length(aco_in),22);
    for j = 1:22
        
        data(:,j) = ((aco_in(:,j+5)-mean(aco_in(:,j+5))))./1E6;
        
    end
    
    timestamp = 1457848722.58 + first_file*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
    
    [psd_val(index+1,:),f] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);
    tlist(index+1) = (index+1)*30;

%     psdplot = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);
%     xlim([5 6000]);
%     ylim([10 110]);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     psd_mat(index) = getframe(gcf);
%     
%     [~,ind0] = min(abs(freq-40));
%     [~,ind1] = min(abs(freq-300));
%     [~,ind2] = min(abs(freq-1500));
%     
%     psdplot_fil1 = psdplot(ind0:ind1);
%     freq_fil1 = freq(ind0:ind1);
%     
%     psdplot_fil2 = psdplot(ind1:ind2);
%     freq_fil2 = freq(ind1:ind2);
%     
%     psd_fit(index+1).file_name = ['ACO',num2str(first_file)];
%     
%     % Fit 40-300Hz
%     figure
%     plot(freq_fil1,10*log10(psdplot_fil1/1E-12),'o','linewidth',1.5)
%     hold on
%     set(gca,'Fontsize',20);
%     title(['PSD Fit; Time = ',data_name,'; Window size = ', num2str(window_size)]);
%     xlabel('Frequency (Hz)');
%     xlim([40 300])
%     ylabel('Power/Freq (dB/Hz)');
%     ylim([40 100]);
%     
%     [p,s] = polyfit(freq_fil1,10*log10(psdplot_fil1/1E-12),1);
%     r2 = 1 - s.normr^2 / norm(10*log10(psdplot_fil1/1E-12)-mean(10*log10(psdplot_fil1/1E-12)))^2;
%     plot(freq_fil1,p(1)*freq_fil1+p(2),'r-','linewidth',1.5)
%     equation = ['y = ',num2str(p(1)),'x + ',num2str(p(2)),';   r^2 = ',num2str(r2)];
%     text(200,80,equation,'FontSize',20);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     saveas(gcf,[pwd ['/band1_psd_fit/T0055853_ACO',num2str(first_file),'.png']]);
%     
%     psd_fit(index+1).band1_slope = p(1);
%     psd_fit(index+1).band1_yint = p(2);
%     psd_fit(index+1).band1_r2 = r2;
%     
%     
%     
%     % Fit 300-1500Hz
%     figure
%     plot(freq_fil2,10*log10(psdplot_fil2/1E-12),'o','linewidth',1.5)
%     hold on
%     set(gca,'Fontsize',20);
%     title(['PSD Fit; Time = ',data_name,'; Window size = ', num2str(window_size)]);
%     xlabel('Frequency (Hz)');
%     xlim([300 1500])
%     ylabel('Power/Freq (dB/Hz)');
%     ylim([20 90]);
%     
%     [p,s] = polyfit(freq_fil2,10*log10(psdplot_fil2/1E-12),1);
%     r2 = 1 - s.normr^2 / norm(10*log10(psdplot_fil2/1E-12)-mean(10*log10(psdplot_fil2/1E-12)))^2;
%     plot(freq_fil2,p(1)*freq_fil2+p(2),'r-')
%     equation = ['y = ',num2str(p(1)),'x + ',num2str(p(2)),';   r^2 = ',num2str(r2)];
%     text(750,60,equation,'FontSize',20);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     saveas(gcf,[pwd ['/band2_psd_fit/T0055853_ACO',num2str(first_file),'.png']]);
%     
%     psd_fit(index+1).band2_slope = p(1);
%     psd_fit(index+1).band2_yint = p(2);
%     psd_fit(index+1).band2_r2 = r2;
%     
% 
end

tlist = [0 tlist];
% v_psd = VideoWriter(['T0055853_PSD_ACO',num2str(first_file),'_ti800s.avi','Uncompressed AVI']);
% v_psd.FrameRate = 3;
% open(v_psd)
% writeVideo(v_psd,psd_mat)
% close(v_psd)

%% Normplot dt = 2min

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
data(:,1) = third_oband(:,4);
data(:,2) = third_oband(:,6);
data(:,3) = third_oband(:,7);
data(:,4) = third_oband(:,8);
data(:,5) = third_oband(:,10);
data(:,6) = third_oband(:,13);
data(:,7) = third_oband(:,15);
data(:,8) = third_oband(:,17);
data(:,9) = third_oband(:,20);
data(:,10) = third_oband(:,22);
data(:,11) = third_oband(:,24);

figure
h = normplot(data);
set(h,'linewidth',3);
set(gca,'Fontsize',20);
xlabel('dB');
legend({'40Hz','63Hz','80Hz','100Hz','160Hz','315Hz','500Hz','800Hz','1600Hz','2500Hz','4000Hz'});

mean_list = [];
std_list = [];

for i = 1:size(data,2)
    mean_list(i) = mean(data(:,i));
    std_list(i) = std(data(:,i));
end

figure
plot(tlist,data(:,1))
hold on
plot(tlist,data(:,4))
hold on
plot(tlist,data(:,7))
hold on
plot(tlist,data(:,8))
hold on
plot(tlist,data(:,10))
set(gca,'Fontsize',20);
xlabel('Time (s)');
ylabel('dB');
legend({'40Hz','100Hz','500Hz','800Hz','2500Hz'});

%% Freq Correlation

[~,ind1] = min(abs(f-40));
[~,ind2] = min(abs(f-60));
[~,ind3] = min(abs(f-80));
[~,ind4] = min(abs(f-100));
[~,ind5] = min(abs(f-300));
[~,ind6] = min(abs(f-500));
[~,ind7] = min(abs(f-800));
[~,ind8] = min(abs(f-1000));
[~,ind9] = min(abs(f-1500));
[~,ind10] = min(abs(f-2500));
[~,ind11] = min(abs(f-4000));

indlist = [ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,ind10,ind11];

r = zeros(length(indlist),size(psd_val,2));
for i = 1:length(indlist)
    for j = 1:size(psd_val,2)
        cov_mat = cov(psd_val(:,indlist(i)),psd_val(:,j));
        r(i,j) = cov_mat(1,2)/(sqrt(cov_mat(1,1))*sqrt(cov_mat(2,2)));
    end
end
 
freqlist = [40 60 80 100 300 500 800 1000 1500 2500 4000];
figure
for ind = 1:11
h = plot(f,r(ind,:));
grid on
set(h,'linewidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Correlation Coefficient');
title([num2str(freqlist(ind)),' Hz']);
set(gca,'Fontsize',20);
saveas(gcf,[pwd ['/f_correlation/',num2str(freqlist(ind)),'Hz.fig']]);
end

