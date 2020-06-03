FS = 12000;     
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 11;
window_size = 4096;
window = hanning(window_size);
overlap = 0.5;
NFFT = 2*window_size;

set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = length (directory);
first_file = 1000;
last_file = 17000;
psd_val = [];

% Read DATA
aco_in = zeros (NUM_SAMPLES, 32);

% Start looping over ACO*.DAT files
for i = first_file:last_file

    disp([num2str(i),' / ',num2str(last_file)])

    filename = [prefix directory(i).name];
    fid = fopen (filename, 'r', 'ieee-le');

    if (fid <= 0)
        continue;
    end

    % Read the single precision float acoustic data samples (in uPa)
    for j = 1:NUM_CHANNELS
        aco_in(:,j) = fread (fid, NUM_SAMPLES, 'float32');
    end

    fclose (fid);

    % Nomalized to zero mean; take middle 22 channels
    aco_norm = zeros(length(aco_in),22);
    for n = 1:22
        aco_norm(:,n) = ((aco_in(:,n+5)-mean(aco_in(:,n+5))))./10^6;
    end

    timestamp = 1457848722.58 + i*2;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);

   [psd_val(i-first_file+1,:),f] = psd(aco_norm,window_size,window,overlap,NFFT,FS,chn,data_name);
%     psd(aco_norm,window_size,window,overlap,NFFT,FS,chn,data_name);
%     xlim([5 6000]);
%     ylim([10 110]);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     psd_mat(i-first_file+1) = getframe(gcf);

end

% v_psd = VideoWriter('ICEX16_PSD_ACO1000_2000_ti2s.avi','Uncompressed AVI');
% v_psd.FrameRate = 3;
% open(v_psd)
% writeVideo(v_psd,psd_mat)
% close(v_psd)

%% Sound Variability Normplot dt = 2s

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
data(:,2) = third_oband(:,7);
data(:,3) = third_oband(:,10);
data(:,4) = third_oband(:,13);
data(:,5) = third_oband(:,15);
data(:,6) = third_oband(:,17);
data(:,7) = third_oband(:,19);
data(:,8) = third_oband(:,22);
data(:,9) = third_oband(:,24);


figure
normplot(data)
legend({'40Hz','80Hz','160Hz','315Hz','500Hz','800Hz','1250Hz','2500Hz','4000Hz'});
        