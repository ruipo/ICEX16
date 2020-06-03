

% INPUTS
%
% data		time series (time x depth)
% p		element positions matrix,  [x_e1,y_e1,z_e1;x_e2...]
% FS		sample frequency (Hz)
% elev		desired look angles in elevation in degrees
% az    desired look angles in azimuth in degrees
% c		sound speed
% f_range     frequency range of interest
% window    frequency window weighting
% overlap    percent overlap [0 1]
% weighting     spatial array weighting
% NFFT
%
% OUTPUTS
%
% MPDR beamform_output	beamformed result (time x elev x az x freq)

function[beamform_output,t,t_end,Sx_mat] = beamform_MPDR_3D(data,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting)

% Define variables
N = size(p,1);
beam_elev = (90-elev).*(pi/180);
beam_az = az.*(pi/180);

win_len = length(window);
t_end = size(data,1)/FS;
df = FS/NFFT;
flist = (-FS/2:df:FS/2-df)';

for m = 2:size(data,2)
    window = [window,window(:,1)];
end

% Format data
data_len = 2^nextpow2(size(data,1));
zero_len = data_len - size(data,1);
zero_mat = zeros(zero_len,size(data,2));
data = [data;zero_mat];

window_start = round(win_len-win_len*overlap);
num_window = round(size(data,1)/window_start)-1;
beamform_output = zeros(num_window,length(beam_elev),length(beam_az),length(f_range));
t = zeros(num_window,1);

% FFT Data
ts_f_mat = zeros(NFFT,N,num_window);
for l = 1:num_window
    ts_f_mat(:,:,l) = fftshift(fft(window.*data(l*window_start-window_start+1:l*window_start-window_start+win_len,:),NFFT));
    t(l) = ((l+1)*window_start-window_start+1)/FS;
end

% Calculate Sn_mat for plane wave noise
Sx_mat = zeros(N,N,length(f_range));
for ii = 1:length(f_range) 
    %ii
    f_val = flist-f_range(ii);
    [~,ind] = min(abs(f_val));
    ts_ff = squeeze(ts_f_mat(ind,:,:));
    Sx_mat(:,:,ii) = (ts_ff * ts_ff')/num_window;
end
        

% for ii = 1:length(f_range) 
%     ii
%     Sx_t = 0;
%     for l = 1:num_window
%         ts_f = fftshift(fft(window.*data(l*window_start-window_start+1:l*window_start-window_start+win_len,:),NFFT));
%         t(l) = ((l+1)*window_start-window_start+1)/FS;
% 
%         f_val = flist-f_range(ii);
%         [~,ind] = min(abs(f_val));
% 
%         ts_f = ts_f(ind,:).';
%         Sx_t = Sx_t + ts_f*ts_f';
%     end
%      
%     Sx_mat(:,:,ii) = Sx_t/num_window;
% end
 
% Start beamforming
%for l = 1:num_window
    %disp([num2str(l),' / ' num2str(num_window)])
    %ts_f = fftshift(fft(window.*data(l*window_start-window_start+1:l*window_start-window_start+win_len,:),NFFT));
    %t(l) = ((l+1)*window_start-window_start+1)/FS;
    
    f_val_start = flist-f_range(1);
    [~,ind_start] = min(abs(f_val_start));

    f_val_end = flist-f_range(end);
    [~,ind_end] = min(abs(f_val_end));
    
    ts_f_mat = ts_f_mat(ind_start:ind_end,:,:);
    
    % calculate k
    k = 2*pi*flist(ind_start:ind_end)./c;

    % linear window
    if strcmp(weighting,'uniform')
        win = ones(1,N)./N;
    end

    % build steering vectors
    for j = 1:length(beam_az)
        for mm = 1:length(beam_elev)
            %mm
            % form MVDR steering vector
            steer = exp(1i * k *(sin(beam_elev(mm))*cos(beam_az(j))*p(:,1)'+sin(beam_elev(mm))*sin(beam_az(j))*p(:,2)'+cos(beam_elev(mm))*p(:,3)'));
            steer_MPDR = zeros(size(steer));
            
            for ii = 1:length(f_range)
                steer_f = steer(ii,:).';
                steer_MPDR(ii,:) = (Sx_mat(:,:,ii)\steer_f) / (steer_f'*(Sx_mat(:,:,ii)\steer_f));   
            end

            % apply weighting
            steer_MPDR = steer_MPDR.*(ones(size(k,1),1)*win);
            
            for l = 1:num_window
                % beamform
                b_elem = sum((conj(steer_MPDR).*ts_f_mat(:,:,l)).');

                beamform_output(l,mm,j,:) = b_elem./(N*NFFT);
            end
        end
    end
%end

t = t - t(1);
t_end = t_end - t(1);
end

