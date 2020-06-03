FS = 12000;
NUM_SAMPLES = FS*2;     
NUM_CHANNELS = 32;
elev = -90:1:90;
z = -1*[7.875,7.125,6.375,5.625,4.875,4.125,3.375,2.625,1.875,1.125,0.375,-0.375,-1.125,-1.875,-2.625,-3.375,-4.125,-4.875,-5.625,-6.375,-7.125,-7.875];
window = hamming(4096);
overlap = 0.5;
NFFT = 2*length(window);
taper = ones(length(z),1)./22;
f_range = [800 900];
c_0 = 1435;
adaptive = 0;
btr = 1;
lofar = 0;
fraz = 0;

% Set Path to DATA
prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-13/DURIP/DURIP_20160313T055853/';
%prefix = '/Volumes/icex6/ICEX_UNCLASS/ICEX16/macrura/2016-03-14_andbefore/DURIP/DURIP_20160314T002324/';
%prefix = '/Users/Rui/Documents/Graduate/Research/ICEX16/ICEX_data/test/';

directory = dir([prefix 'ACO0000*.DAT']);
num_files = 15;
last_file = 1060;

for index = 0:(last_file-1045)/num_files-1
  
    disp([num2str(index),' / ', num2str((last_file-1045)/num_files-1)])
    
    first_file = 1045 + index*num_files;

    % Read DATA
    aco_in = zeros (NUM_SAMPLES * num_files, 32);

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
    aco_norm = zeros(length(aco_in),22);
    for i = 1:22
        
        aco_norm(:,i) = ((aco_in(:,i+5)-mean(aco_in(:,i+5))))./10^6;
        
    end
    
    clear aco_in
    
    timestamp = 1457848722.58 + first_file*2;
    %timestamp = 1457914993.67 + first_file*2;
    data_name = ['ACO',num2str(first_file)];
    
    [beamform_elev,beamform_elev_allf,f,t] = vert_array_beamform(aco_norm,elev,z,window,overlap,NFFT,FS,taper,f_range,c_0,adaptive,btr,lofar,fraz,data_name);
    saveas(gcf,[pwd ['/btr/ACO',num2str(first_file),'.fig']]);
    
end





