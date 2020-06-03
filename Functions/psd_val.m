
% Returns the power spectral density value of the input data at a specific frequency or averaged over a frequecy band.
% aco_in must be [# of data points] by [# of channels].
% Can input just one channel or multiple channels in a list.
% Can input just one frequency or a frequency band as [start freq, end freq].

function psd_val = psd_val(data,window,NFFT,f_range,FS,channel,data_name)

ts = data(:,channel);
[s,f,t] = spectrogram(ts,window,[],NFFT,FS);

if length(f_range) > 1
    [~,ind_start] = min(abs(f-f_range(1)));
    [~,ind_end] = min(abs(f-f_range(end)));

    psd_val = mean((1/(FS*NFFT)) * abs(s(ind_start:ind_end,:)).^2,1);
end

if length(f_range) == 1
    [~,ind] = min(abs(f-f_range));
    psd_val = (1/(FS*NFFT)) * abs(s(ind,:)).^2;
end

% 
% figure
% plot(t,10*log10(psd_val/1E-6^2));
% set(gca,'Fontsize',20);
% title(['PSD vs. Time; ',num2str(f_range(1)),' - ',num2str(f_range(end)),' Hz Band; Timestamp = ',data_name]);
% xlabel('Time (s)');
% ylabel('Power/Frequency (dB/Hz)');

end
