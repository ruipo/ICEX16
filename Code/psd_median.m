clear
FS = 12000; 
window_size = 2048;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 1;
data_name = 'Hour 1';

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
directory = dir([prefix 'ACO0000*.DAT']);

for iter = 1
    iter
first_file = 2000 + (1800*iter - 1800);
last_file = first_file+1;

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

[psd_med(iter,:),freq] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);

 semilogx(freq,10*log10(psd_med(iter,:)/(1E-6)^2),'linewidth',1.5)
 grid on
 hold on

end

