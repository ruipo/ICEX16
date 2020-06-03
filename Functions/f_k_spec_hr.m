%returns the frequency wavenumber spectrum of the input data
%requires data in timeseries X channel

function [P,freq,k,cov_mat] = f_k_spec_hr(data,window,overlap,NFFT,FS,c,z)
        
        N = size(data,2);
        win_len = NFFT;

        window_start = round(win_len-win_len*overlap);
        num_window = round(size(data,1)/window_start)-2;
        t = zeros(num_window,1);
        
        freq = 0:FS/NFFT:FS/2;
        
        for m = 2:size(data,2)
            window = [window,window(:,1)];
        end

        % FFT Data
        ts_f_mat = zeros(NFFT,N,num_window);
        ts_f_mat_onesided = zeros(NFFT/2+1,N,num_window);
        for l = 1:num_window
            ts_f_mat(:,:,l) = fft(window.*data(l*window_start-window_start+1:l*window_start-window_start+win_len,:),NFFT);
            ts_f_mat_onesided(:,:,l) = (1/(sqrt(FS)*norm(window(:,1),2))).*ts_f_mat(1:NFFT/2+1,:,l);
            ts_f_mat_onesided(2:end-1,:,l) = sqrt(2).*ts_f_mat_onesided(2:end-1,:,l);
            t(l) = ((l+1)*window_start-window_start+1)/FS;
        end

        % Calculate Sn_mat for plane wave noise
        Qx_mat = zeros(N,N,length(freq));
        for ii = 1:length(freq) 
            ts_ff = squeeze(ts_f_mat_onesided(ii,:,:));
            cov_mat(:,:,ii) = (ts_ff * ts_ff')/num_window;
            Qx_mat(:,:,ii) = inv((ts_ff * ts_ff')/num_window);
        end 
        
        alpha = 1/c;
        k1 = freq*alpha;
        k2 = -fliplr(k1);
        k2(end) = [];
        k = [k2 k1];

        Pinv = zeros(length(freq),length(k));
        P = zeros(length(freq),length(k));
        
        for f = 1:length(freq)
            f
            for n = 1:N
                for m = 1:N
                    Pinv(f,:) = Pinv(f,:) + Qx_mat(m,n,f)*exp(1i*2*pi*k*(z(n)-z(m)));
                end
            end
            P(f,:) = 1./Pinv(f,:);
        end
        
        %P = (1/N^2).*P;
                
end
