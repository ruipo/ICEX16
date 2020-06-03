%% Import time and AEL info
path = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/icex_spl_sd/alogs/macrura_13_3_2016_____05_58_43_alvtmp/';

ael_DB_time = importdata([path 'DB_TIME.klog']);
ael_uptime = str2double(ael_DB_time.textdata(:,1));
ael_db_utc = ael_DB_time.data;
ael_db_utc_0 = ael_db_utc(1)-ael_uptime(1);

ael_pitch_t1 = importdata([path 'AEL_PITCH_C.klog']);
ael_pitch_t2 = sprintf('%s*', (ael_pitch_t1.textdata{:}));
ael_pitch = ael_pitch_t1.data; 
ael_pitch_uptime = str2double(regexp(ael_pitch_t2,'[+-]?\d+\.?\d*','match')).';
ael_pitch_utc = ael_pitch_uptime + ael_db_utc_0;
clear('ael_pitch_t1','ael_pitch_t2');

ael_depth_t1 = importdata([path 'NAV_Z.klog']);
ael_depth_t2 = sprintf('%s*', (ael_depth_t1.textdata{:}));
ael_depth = ael_depth_t1.data; %auv z
ael_depth = ael_depth - 38;
ael_depth_uptime = str2double(regexp(ael_depth_t2,'[+-]?\d+\.?\d*','match')).';
ael_depth_utc = ael_depth_uptime + ael_db_utc_0; %auv z utc time
clear('ael_depth_t1','ael_depth_t2');

%% Import data

FS = 12000; 
window_size = 2048;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 16;

starttime_str = '2016-03-13 05:58:53';
starttime = datetime(starttime_str);
starttime_utc = posixtime(starttime);

% Element Locations
N = 32;
% 0.75m spacing in the middle
N2 = 22;
z2 = zeros(N2,1);
d = 0.75;
for n = 0:N2-1
    z2(n+1) = -(n-(N2-1)/2)*d; 
end

% 1.5m spacing on top
N1 = 5;
z1 = zeros(N1,1);
d = 1.5;
z1(end) = z2(1) + d;

for n = 1:N1-1
    z1(end-n) = z1(end-(n-1)) + d;
end

% 1.5m spacing on bottom
N3 = 5;
z3 = zeros(N3,1);
d = 1.5;
z3(1) = z2(end) - d;

for n = 2:N3
    z3(n) = z3(n-1) - d;
end

z = [z1; z2; z3];
p = [zeros(1,N) ; zeros(1,N) ; z'];
p = p';

elev = -90:1:90;
az = 0;
c = 1435;
windowb = hanning(8192);
NFFTb = 1024;
f_range = [100 1000];
weighting = 'icex_hanning';

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';

directory = dir([prefix 'ACO0000*.DAT']);

num_files = 15;
last_file = 17002;


% Read DATAs
aco_in = zeros (NUM_SAMPLES * num_files, 32);
time_vec = zeros(1,round(last_file/num_files));
psd_mat = zeros(round(last_file/num_files),overlap*NFFT+1);
psd_var_mat = zeros(round(last_file/num_files),overlap*NFFT+1);
beamform_mat = zeros(round(last_file/num_files),length(elev));
beamform_f_mat = zeros(length(elev),NFFTb,round(last_file/num_files));

for index = 0:round(last_file/num_files)-2
    
    disp([num2str(index),' / ', num2str(round(last_file/num_files)-1)])
    
    first_file = 2000 + index*num_files;

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

    % Nomalized to zero mean;
    data = zeros(length(aco_in),32);
    for j = 1:32       
        data(:,j) = ((aco_in(:,j)-mean(aco_in(:,j))))./10^6;
    end
    
    time_vec(index+1) = (first_file-1)*2 + (num_files/2)*2;
    
%     % PSD
%     [~,psd_avg,psd_var,freq] = psd(data,window_size,window,overlap,NFFT,FS,16,'test');
%     
%     psd_mat(index+1,:) = psd_avg;
%     psd_var_mat(index+1,:) = psd_var;
    
    % BEAMFORMING
    t_pitch = (first_file+num_files-1)+starttime_utc;
    pitch = interp1(ael_pitch_utc,ael_pitch,t_pitch);
    p_pitch = p;
    p_pitch(:,1) = p(:,3)*cos(deg2rad(90-pitch));
    p_pitch(:,3) = p(:,3)*sin(deg2rad(90-pitch));

    [beamform_output,t,t_end] = beamform_3D(data,p_pitch,FS,elev,az,c,f_range,NFFTb,windowb,overlap,weighting);
    beamform_f = squeeze(mean(beamform_output,1));
    beamform = squeeze(mean(mean(beamform_output,1),4));
    
    beamform_mat(index+1,:) = beamform;
    beamform_f_mat(:,:,index+1) = beamform_f;

end

%%
psd_mat(end-62:end,:) = [];
beamform_mat(end-62:end,:) = [];
time_vec(end-62:end) = [];
beamform_f_mat(:,:,end-62:end) = [];

timevec_utc = time_vec+starttime_utc;
datetimevec = datetime(timevec_utc, 'convertfrom','posixtime');
depth_interp = interp1(ael_depth_utc,ael_depth,timevec_utc);

%% Plotting - time vs freq @ bf angle

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/icex_spl_sd/03-13-2016/ang_figures/';

fr = linspace(100,1000,NFFTb);
f1 = 800;
f2 = 900;
[~,ind1] = min(abs(fr-f1));
[~,ind2] = min(abs(fr-f2));

bf_all = squeeze(mean(median(beamform_f_mat(:,ind1:ind2,:),3),2));

figure('units','normalized','outerposition',[0 0 1 1]);
s1 = subplot(3,3,[1,2,4,5]);
ylabel('Frequency (Hz)');
xlabel('Time (03-13-2016 UTC)')
datetick('x', 'HH:MM','keeplimits');
ylim([100 1000])
xlim([datenum(datetimevec(1)) datenum(datetimevec(end))])
%Xticks([])
%set(gca,'Xticklabel',[]);
%set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([30,90])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
hold on

s2 = subplot(3,3,[7,8]);
plot(datenum(datetimevec),-depth_interp,'k','linewidth',2)
ylim([35 40])
datetick('x', 'HH:MM','keeplimits');
xlim([datenum(datetimevec(1)) datenum(datetimevec(end))])
set(gca,'Ydir','reverse')
xlabel('Time (03-13-2016 UTC)')
ylabel('Depth (m)')
set(gca,'fontsize',15)
grid on
hold on

s3 = subplot(3,3,[3,6,9]);
plot(10*log10(bf_all/(1E-6)^2),elev,'b','linewidth',2)
ylim([-90 90])
xlim([30 90])
yticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Yticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]) 
grid on
ylabel('Elevation (Degrees)');
xlabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
title('Frequency = 850 Hz')
hold on


for ee = 1:length(elev)
    ee
    figure(1)
    axes(s1);
    h = pcolor(datenum(datetimevec),fr,10*log10(squeeze(beamform_f_mat(ee,:,:))/(1E-6)^2));
    set(h,'Edgecolor','None')
    datetick('x', 'HH:MM','keeplimits');
    ylim([100 1000])
    xlim([datenum(datetimevec(1)) datenum(datetimevec(end))])
    title(['Elevation = ' num2str(elev(ee)) ' Degrees'])

    axes(s3);
    f2 = plot([30 90],ones(2,1)*elev(ee), 'k','linewidth',1.5);
    
    saveas(gcf,[figpath num2str(ee) '.png']);
    pause(0.1)
    
    delete(h);
    delete(f2);
    
end

%% make video - ANG

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'ang.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);

%% plotting - Time vs dB @ frequency
figure

flist = [100 500 1000 1500 2000 3000 6000];
color = ['k','r','m','y','g','b','c'];
shape = ['-','-','-','-','-','-','-'];
for ff = 1:length(flist)
    f = flist(ff);
    [~,ind] = min(abs(f-freq));

    plot(datetimevec, 10*log10(psd_mat(:,ind)/(1E-6)^2),[color(ff) shape(ff)],'linewidth',1)
    hold on
    title('Noise Level vs. Depth');
    ylabel('NL (dB re 1\muPa^2/Hz)')
    xlabel('Time (UTC)')
    grid on
    set(gca,'fontsize',20);
end

legend('100Hz', '500Hz','1000Hz','1500Hz','2000Hz','3000Hz','6000Hz');

%% Plotting - NL

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/icex_spl_sd/03-13-2016/figures/';

figure('units','normalized','outerposition',[0 0 1 1]);
s1 = subplot(3,3,[1,4,7]);
ylim([0 250])
set(gca,'Ydir','reverse')
ylabel('Depth (m)')
set(gca,'fontsize',15)
set(gca,'Xticklabel',[])
yticks([0 25 50 75 100 125 150 175 200 225 250])
set(gca,'Yticklabel',[0 25 50 75 100 125 150 175 200 225 250]);
grid on
hold on

s2 = subplot(3,3,[2,3,5,6]);
h = pcolor(freq,datenum(datetimevec),10*log10(psd_mat/(1E-6)^2));
xlabel('Frequency (Hz)');
ylabel('Time (03-13-2016 UTC)')
datetick('y', 'HH:MM','keeplimits');
xlim([40 6000])
set(h,'Edgecolor','None')
set(gca,'XScale','log')
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
xticks([50 100 200 400 800 1600 3200 6000])
set(gca,'Xticklabel',[50 100 200 400 800 1600 3200 6000]) 
caxis([30,90])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
hold on

s3 = subplot(3,3,[8,9]);
xlim([40 6000])
ylim([30 90])
xticks([50 100 200 400 800 1600 3200 6000])
set(gca,'Xticklabel',[50 100 200 400 800 1600 3200 6000]) 
grid on
set(gca,'XScale','log')
xlabel('Frequency (Hz)');
ylabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
hold on


for tt = 1:length(datetimevec)
    tt
    figure(1)
    axes(s1);
    f3 = plot([0 1],ones(2,1)*-depth_interp(tt),'k','linewidth',2);
    f4 = plot(0.5,-depth_interp(tt),'p','Markersize',15,'MarkerEdgeColor','black','MarkerFaceColor','Yellow');
    title(['Depth = ' num2str(round(-depth_interp(tt),2)) ' m'])

    axes(s2);
    f1 = plot([40 6000],ones(2,1)*datenum(datetimevec(tt)),'k','linewidth',2);

    axes(s3);
    f2 = plot(freq,10*log10(psd_mat(tt,:)/(1E-6)^2),'b','linewidth',2);
    
    saveas(gcf,[figpath num2str(tt) '.png']);
    pause(0.1)

    delete(f1);
    delete(f2);
    delete(f3);
    delete(f4);
end


%% make video - NL

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'nl.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);
  
%% Plotting - depth vs bf max
figure
hold on
title('Max NL Elevation Angle vs. Depth (Frequency = 850Hz)');
xlabel('Elevation (Degrees)')
ylabel('Time (UTC)')
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
set(gca,'Ydir','reverse')
grid on
set(gca,'fontsize',20);

for ii = 1:size(beamform_mat,1)
    [~, inds] = max(10*log10(abs(beamform_mat(ii,:))/1E-12));
    
    plot(elev(inds),datetimevec(ii),'bo','linewidth',1.5);
end
    
%% Plotting - depth vs bf

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/ICEX16/icex_spl_sd/03-13-2016/bf_figures/';

figure('units','normalized','outerposition',[0 0 1 1]);
s1 = subplot(3,3,[1,4,7]);
ylim([0 250])
set(gca,'Ydir','reverse')
ylabel('Depth (m)')
set(gca,'fontsize',15)
set(gca,'Xticklabel',[])
yticks([0 25 50 75 100 125 150 175 200 225 250])
set(gca,'Yticklabel',[0 25 50 75 100 125 150 175 200 225 250]);
grid on
hold on

s2 = subplot(3,3,[2,3,5,6]);
h = pcolor(elev,datenum(datetimevec),10*log10(abs(beamform_mat)/1E-12));
xlabel('Elevation (Degrees)');
ylabel('Time (03-13-2016 UTC)')
datetick('y', 'HH:MM','keeplimits');
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
set(h,'Edgecolor','None')
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([30,90])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
hold on
title('Frequency = 850Hz')

s3 = subplot(3,3,[8,9]);
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
ylim([30,90])
grid on
xlabel('Elevation (Degrees)');
ylabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
hold on


for tt = 1:length(datetimevec)
    tt
    figure(1)
    axes(s1);
    f3 = plot([0 1],ones(2,1)*-depth_interp(tt),'k','linewidth',2);
    f4 = plot(0.5,-depth_interp(tt),'p','Markersize',15,'MarkerEdgeColor','black','MarkerFaceColor','Yellow');
    title(['Depth = ' num2str(round(-depth_interp(tt),2)) ' m'])

    axes(s2);
    f1 = plot([-90 90],ones(2,1)*datenum(datetimevec(tt)),'k','linewidth',2);

    axes(s3);
    f2 = plot(elev,10*log10(abs(beamform_mat(tt,:))/1E-12),'b','linewidth',2);
    
    saveas(gcf,[figpath num2str(tt) '.png']);
    pause(0.1)

    delete(f1);
    delete(f2);
    delete(f3);
    delete(f4);
end

%% make video - BF

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'bf.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);




