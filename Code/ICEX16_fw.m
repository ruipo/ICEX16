
FS = 12000; 
window_size = 4096;
window = hanning(window_size);
NFFT = window_size;
overlap = 0.5;

NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
chn = 1;
data_name = 'Hour 1';

% Set Path to DATA
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
directory = dir([prefix 'ACO0000*.DAT']);

%for p = 6:8
%first_file = 2000 + (1800*p - 1800);
%last_file = first_file+180;

first_file = 5000;
last_file = first_file+360;

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
data = aco_in./1E6;

c = 1435;
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


[P,freq,k,cov_mat] = f_k_spec_hr(data,window,overlap,NFFT,FS,c,z);

%P_mat(:,:,p) = P;
%end

%%

%P = mean(P_mat,3);
figure
fig = pcolor(k,freq,10*log10(abs(P)/1E-12));
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Wavenumber')
ylabel('Frequency (Hz)')
colorbar;
colormap 'jet';
xlim([-0.7 0.7]);
ylim([0 1000]);
caxis([30 75]);

hold on
f1list = fliplr(freq);
f1list(end) = [];
flist = [f1list,freq];
plot(k,flist)
%%
finterest = 1:3:400;

for i = 1:length(finterest)
    
flist = freq - finterest(i);
[~,loc] = min(abs(flist));
kinterest = freq(loc)/c;
dk = k(2)-k(1);

klist = k-(-0.7);
[~,loc1] = min(abs(klist));

klist = k-(-kinterest);
[~,loc2] = min(abs(klist));

klist = k-kinterest;
[~,loc3] = min(abs(klist));

klist = k-(0.7);
[~,loc4] = min(abs(klist));

P_outside(i) = 10*log10(mean(mean(abs(P(loc,loc1:loc2))) + mean(abs(P(loc,loc3:loc4))))/1E-12);
P_cone(i) = 10*log10(mean(abs(P(loc,loc2:loc3)))/1E-12);
end

figure
plot(finterest,P_cone)
hold on
plot(finterest,P_outside)
hold on
plot(finterest,P_cone-P_outside)





