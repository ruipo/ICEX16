%% Set parameters

% data read-in parameters
FS = 12000; % sampling frequency
NUM_SAMPLES = FS*2; % number of samples per file
NUM_CHANNELS = 32; % number of channels

% beamform parameters
array_orientation = 'h'; %'v' for vertical, 'h' for horizontal
c = 1435; % assumed sound speed
window = hanning(8192); % fft window
NFFT = 512; % number of fft frequency bins
f_range = [800 900]; % frequency range to beamform over
overlap = 0.5; % overlap between successive samples (0 to 1)
weighting = 'icex_hanning'; % beamform window weighting


%% Load Data
% Set Path to DATA
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/'; %vertical test
 prefix = '/Users/Rui/Documents/Graduate/Research/ICEX_SIMI/engin_test/DURIP_12_12_2019T193012/'; %horizontal test 
directory = dir([prefix 'ACO0000*.DAT']);

first_file = 1; % first file to read in
last_file = first_file+400; % last file to read in (1800 files = 1hr)
 
% Read DATA
aco_in = zeros(NUM_SAMPLES * (last_file-first_file), 32);
  
% Start looping over ACO*.DAT files
counter=0;
disp('Reading in data ...');
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

data = (aco_in-mean(aco_in,1))./10^6; % de-mean and convert data to uPa
time = 0:1/FS:size(data,1)/FS-1/FS; % create time axis


%% Beamforming

[p,elev,az] = getICEXArray(array_orientation); % get array element locations and beamforming elev and az angles

[beamform_output,t,t_end] = beamform_3D(data,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting); % beamform_ouput [time,elev,az,frequency]

[~,ind] = min(abs(t-t_end));
t = t(1:ind-1);
beamform_output = beamform_output(1:ind-1,:,:,:);
beamform_output_db = squeeze(10*log10(abs(beamform_output)/1E-12)); % beamform_ouput_db [time,ang(elev or az),frequency] units: db re 1uPa

beamform_ang_db = mean(beamform_output_db,3); % beamform_ang_db [time,ang(elev or az)] units: db re 1uPa
beamform_avgt_db = mean(beamform_ang_db,1); % beamform_avgt_db [ang(elev or az)] units: db re 1uPa

%% Plot Figure
if array_orientation == 'v'
    angstr = 'Elevation (Degrees)';
    ang = elev;
else
    angstr = 'Azimuth (Degrees)';
    ang = az;
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
h1 = pcolor(ang,t,beamform_ang_db);
set(h1,'Edgecolor','None')
colorbar
ylabel('Time (s)')
xlabel(angstr)
title(['Frequency = ' num2str(mean(f_range)) 'Hz'])
set(gca,'fontsize',25);
colormap jet

subplot(1,2,2)
plot(beamform_avgt_db,ang,'linewidth',2);
xlabel('dB re 1\muPa');
ylabel(angstr);
ylim([min(ang),max(ang)])
title('Time-averaged Output')
grid on
set(gca,'fontsize',25);

% save figure
savename = ['beamform_' num2str(first_file) '_' num2str(last_file) '.m'];
saveas(gcf,[savename '.fig'])
saveas(gcf,[savename '.png'])