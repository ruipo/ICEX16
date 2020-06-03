clear
FS = 12000; 
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 16;
 
% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';

bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',20,'CutoffFrequency2',1500,'SampleRate',12000);
 
directory = dir([prefix 'ACO0000*.DAT']);

for n = 1
    n
%first_file = randsample(14000,1) + 2000;
first_file = 2000 + (1800*n - 1800);
last_file = first_file + 1800;

startlist(n) = first_file;
 
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

aco_in = filtfilt(bandfilt,aco_in(:,chn));


 
    % Nomalized to zero mean; take selected channel
    data = (aco_in-mean(aco_in))./10^6;
    
    %figure
    %normplot(data)

    dt = 2048;
    overlap = 0.5;

    window_start = round(dt-dt*overlap);
    num_window = round(size(data,1)/window_start)-2;

    t = zeros(num_window,1);
    ts = zeros(dt,num_window);
    for l = 1:num_window
        ts(:,l) = data(l*window_start-window_start+1:l*window_start-window_start+dt);
        t(l) = ((l+1)*window_start-window_start+1)/FS;
    end
    meanlist = mean(ts);
    varlist = var(ts);
    kurtlist = kurtosis(ts);
    skewlist = skewness(ts);
    
end

%%
FS2 = 1/(t(1));
band_100hz = third_oband(19,:,1);
dt2 = 128;
overlap = 0.5;

window_start = round(dt2-dt2*overlap);
num_window = round(length(band_20hz)/window_start)-2;
t2 = zeros(num_window,1);
ts = zeros(dt2,num_window);
    for l = 1:num_window
        ts(:,l) = band_20hz(l*window_start-window_start+1:l*window_start-window_start+dt2);
        t2(l) = ((l+1)*window_start-window_start+1)/FS2;
    end

kurtlist_100hz = kurtosis(10*log10(ts./1E-12));

figure
plot(t2/60,kurtlist_100hz);
xlabel('Time (min)')
ylabel('Kurtosis');
set(gca,'fontsize',25);
grid on

%% Plotting

figure
plot(t/60,gradient(meanlist));
xlabel('Time (min)')
ylabel('Mean');
set(gca,'fontsize',25);
grid on


figure
plot(t/60,varlist);
xlabel('Time (min)')
ylabel('variance');
set(gca,'fontsize',25);
grid on

figure
plot(t/60,kurtlist);
hold on
plot(t/60,3*ones(length(t),1),'r','linewidth',2)
xlabel('Time (min)')
ylabel('Kurtosis');
set(gca,'fontsize',25);
grid on

figure
plot(t/60,skewlist);
xlabel('Time (min)')
ylabel('Kurtosis');
set(gca,'fontsize',25);
grid on
