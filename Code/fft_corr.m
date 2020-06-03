clear
FS = 12000; 
NFFT = 2048;
df = FS/(2*NFFT);
f = -FS/2:df:FS/2-df;
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 16;

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';

bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',40,'CutoffFrequency2',1500,'SampleRate',12000);
 
directory = dir([prefix 'ACO0000*.DAT']);

filelen = [12005,12050,12200,12250,12300,12350,12400,12450,12500,12550,12600,12650,12700,12750,12800,12850,12900,12950,13000];
fcorr = zeros(length(filelen),1);
for iter = 1;%1:length(filelen)
    iter
first_file = 12000;
last_file = filelen(iter);

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
 
% Nomalized to zero mean; take selected channel
 
data = ((aco_in(:,chn)-mean(aco_in(:,chn))))./10^6;

data = filtfilt(bandfilt,data);
time = (1/(FS))*(0:length(data)-1/FS);
 
timestamp = 1457848722.58 + first_file*2;
data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
 
dt = NFFT;
overlap = 0.5;
 
window_start = round(dt-dt*overlap);
num_window = round(size(data,1)/window_start)-2;
ts = zeros(dt,num_window);

for l = 1:num_window
    ts(:,l) = data(l*window_start-window_start+1:l*window_start-window_start+dt);
end

ts_f = fft(ts,2*NFFT);
ts_f = fftshift(ts_f,1);
Sx_sample = abs((ts_f * ts_f'));
PSD = diag(Sx_sample);
multiplier = mean(abs(data).^2)/(sum(PSD)*df);
Sx_sample = multiplier*Sx_sample;
PSD = multiplier*PSD;

plot(f,PSD)

ii = find((f>=-1500)&(f<=1500));
iPSD = sqrt(PSD.^(-1));
normSx = diag(iPSD(ii))*Sx_sample(ii,ii)*diag(iPSD(ii));
fcorr(iter)= (sum(sum(normSx))-length(ii))/(length(ii).^2-length(ii));
    
end