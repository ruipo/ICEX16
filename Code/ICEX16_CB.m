
%% Element Locations
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

%% Load Data
FS = 12000; 
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
beamform_output_db_all = [];
t_all = [];

 
% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/engin_test/DURIP/DURIP_20290817T014032/';
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
 
directory = dir([prefix 'ACO0000*.DAT']);

first_file_list = [2001];

for s = 1
%first_file = 2000 + (1800*s - 1800);
%last_file = first_file + 1800;

first_file = first_file_list(s);
last_file = first_file + 30;
 
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

data = aco_in;
%data = [aco_in(1:10000000,:);aco_in(14000000:24000000,:);aco_in(28000000:38500000,:)];

% convert to Pa
%aco_in = aco_in
% % Nomalized to zero mean; take middle 22 channels
% aco_norm = zeros(length(aco_in),32);
% for i = 1:N
%     aco_norm(:,i) = ((aco_in(:,i)-mean(aco_in(:,i))))./10^6;
% end
% 
% clear aco_in

timestamp = 1457848722.58 + first_file*2;
data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);

%%% MPDR Beamforming

elev = -90:1:90;
%elev = fliplr(elev);
az = 0;
c = 1435;
window = hanning(8192);
NFFT = 512;
%df = FS/NFFT;
f_range = [800 900];
overlap = 0.5;
weighting = 'icex_hanning';

[beamform_output,t,t_end] = beamform_3D(data./10^6,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting);

[~,ind] = min(abs(t-t_end));
t = t(1:ind-1);
beamform_output = beamform_output(1:ind-1,:,:,:);
beamform_output_db = squeeze(10*log10(abs(beamform_output)/1E-12));

t_all = [t_all;(t+t_end*(s-1))];
beamform_output_db_all = [beamform_output_db_all;beamform_output_db];

end

%% Plotting
f = linspace(f_range(1),f_range(2),NFFT);
[~,ind1] = min(abs(f - 800));
[~,ind2] = min(abs(f - 900));

%[~,ind] = min(abs(t-t_end));
beamform_output_t = squeeze(mean(beamform_output_db(1:ind-1,:,:),1)); 
%beamform_output_tmed = squeeze(median(beamform_output_db(1:ind-1,:,:),1));
beamform_800_900hz = squeeze(mean(beamform_output_t(:,ind1:ind2),2));
%beamform_800_900hz_med = squeeze(mean(beamform_output_tmed,2));

figure(2)
hold on
plot(beamform_800_900hz,elev,'linewidth',3)
%hold on
%plot(beamform_800_900hz_med,elev,'linewidth',2)
set(gca,'Fontsize',30);
xlabel('Power (dB re 1\muPa^2)')
ylabel('Elevation (Degrees)')
ylim([-90 90]);
title(['Time = ',data_name])
grid on

%%

beamform_output_f = squeeze(mean(beamform_output_db(:,:,ind1:ind2),3)).'; 

% for n = 1:size(beamform_output_f,1)
%     beamform_output_f(n,1:ind) = (beamform_output_f(n,1:ind)-mean(beamform_output_f(n,1:ind)))/std(beamform_output_f(n,1:ind));
% end

figure
fig = pcolor(t(1:ind-1)./60,elev,beamform_output_f);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (min)')
xlim([0 t_end./60]);
ylabel('Elevation (Degrees)')
title(['Time = ',data_name,'; f = 850Hz'])
colorbar;
colormap 'jet';

%%

figure
fig = pcolor(f_range,elev,beamform_output_t);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Freqnuency (Hz)')
ylabel('Elevation (Degrees)')
title(['Time = ',data_name,'; f = 850Hz'])
colorbar;
colormap 'jet';

%% Plotting
f = linspace(f_range(1),f_range(2),NFFT);
for i = 76;%1:length(elev)
beamform_output_elev = squeeze(beamform_output_db_all(:,i,:)).';

figure
%fig = pcolor(t(1:ind-1)./60,f,beamform_output_elev);
fig = pcolor(t_all./3600,f,beamform_output_elev);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (Hrs)')
xlim([0 t_all(end)/3600]);
ylabel('Frequency (Hz)')
title(['Elevation = ',num2str(elev(i))])
colorbar;
colormap 'jet';
end



