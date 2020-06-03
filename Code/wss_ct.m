clear
FS = 12000; 
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 16;
 
% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';

bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',20,'CutoffFrequency2',1500,'SampleRate',12000);
 
directory = dir([prefix 'ACO0000*.DAT']);

filelen = [5,15,30,60,120,240,480,720,960,1200,1440,1680,1800];
%filelen = [10 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000];

for n = 1:8
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

NFFT = 4*2048;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;

for iter = 1:length(filelen)
    clear data
    iter
 
    % Nomalized to zero mean; take selected channel
    data = ((aco_in(1:NUM_SAMPLES*filelen(iter))-mean(aco_in(1:NUM_SAMPLES*filelen(iter)))))./10^6;
    %time = (1/(FS))*(0:length(data)-1/FS);

    %Rx = abs(autocorr(data,length(data)-1));
    %lagtime(iter) = time(find(Rx >= 0.1,1,'last'));

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
        %t(l) = ((l+1)*window_start-window_start+1)/FS;
    end
    meanlist = mean(ts);
    varlist = var(ts);

    ts_f = fft(ts,NFFT);
    ts_f = fftshift(ts_f,1);
    Sx_sample = abs((ts_f * ts_f'));
    PSD = diag(Sx_sample);
    %multiplier = mean(abs(data).^2)/(sum(PSD)*df);
    multiplier = (1/(FS*NFFT));
    Sx_sample = multiplier*Sx_sample;
    PSD = multiplier*PSD;

    ii = find((f>=-1500)&(f<=1500));
    iPSD = sqrt(PSD.^(-1));
    normSx = diag(iPSD(ii))*Sx_sample(ii,ii)*diag(iPSD(ii));
    fcorr(iter,n)= (sum(sum(normSx))-length(ii))/(length(ii).^2-length(ii));
    
    %figure
    %plot(t,meanlist);

    %figure
    %plot(t,gradient(meanlist));

    meangradlist(iter,n) = median(abs(gradient(meanlist)));
    vargradlist(iter,n) = median(abs(gradient(varlist)));
end
end

%% Plotting

figure
for i = 1:size(meangradlist,2)
plot(filelen*2/60,abs(meangradlist(:,i)),'linewidth',1.5)
hold on
end
xlabel('Time (min)')
ylabel('Median Change in Mean');
set(gca,'fontsize',25);
grid on

figure
for i = 1:size(vargradlist,2)
plot(filelen*2/60,abs(vargradlist(:,i)),'linewidth',1.5)
hold on
end
xlabel('Time (min)')
ylabel('Median Change in Variance');
set(gca,'fontsize',25);
grid on


figure
for i = 1:size(fcorr,2)
plot(filelen*2/60,fcorr(:,i),'linewidth',1.5)
hold on
end
xlabel('Time (min)')
ylabel('Correlation Between Difference Frequency Bins');
set(gca,'fontsize',25);
grid on

