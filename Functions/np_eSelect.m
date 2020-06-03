
% Select for event events with the non-parametric method by searching for
% improbable clusters of high amplitude peaks
% set thresholds w, r, and pthres.

function[peak_event_mat_f,loc_event_mat_f,num_event] = np_eSelect(data,FS,w,r,pthres,view_chn,array_spacing,c_0,data_name)

num_chn = size(data,2);
num_event = zeros(num_chn,1);
loc_event_mat = nan(num_chn,5000);
peak_event_mat = nan(num_chn,5000);

% Event Selection
for chn  = 1:num_chn

    disp([num2str(chn), ' / ',num2str(num_chn)])

    clear data_peaks
    clear loc
    clear width
    clear data_abs
    
    data_abs = abs(data(:,chn));
    loc_event = [];
    peak_event = [];

    [data_peaks,loc] = findpeaks(data_abs);
    [pd, I] = ecdf(data_peaks);
    pdlist = find(abs(pd-(1-r)) < 0.01);
    l = I(pdlist(1));
    
    count = 0;

    for j = 1:length(data_peaks)/w+1
        
        if j*w <= length(data_peaks)
            data_seg = data_peaks(j*w-w+1:j*w);
        end
        
        if j*w > length(data_peaks)
            data_seg = data_peaks(j*w-w+1:end);
        end
        
        n = sum(data_seg>l);

        p = 1 - binocdf(n,length(data_seg),r);

        if p < pthres
            count = count+1;
            
            seg_loc = find(data_seg>l);
            seg_val = data_seg(seg_loc);
            
            for location = 1:length(seg_loc)
                loc_event = [loc_event,loc(j*w-w+seg_loc(location))/FS];
            end
            
            peak_event = [peak_event;seg_val];
        end
        
    end

    loc_event = sort(loc_event);
    num_event(chn) = length(loc_event);
    peak_event_mat(chn,1:length(peak_event)) = peak_event;
    loc_event_mat(chn,1:length(loc_event)) = loc_event;
end

% Retain selection only if its found in atleast half of the channels
loc_event_mat_f = nan(num_chn,5000);
peak_event_mat_f = nan(num_chn,5000);
for channel = 1:num_chn
    clear lia
    clear test
    clear test_peak
    test = loc_event_mat(channel,:);
    test = test(~isnan(test));
    lia_mat = nan(num_chn,length(test));

    for chn = 1:num_chn
        lia_mat(chn,:) = ismembertol(test,loc_event_mat(chn,:),(array_spacing/(c_0))*abs(chn-channel),'DataScale',1);
    end

    lia = sum(lia_mat,1);
    lia(lia <= 3*num_chn/4) = 0;
    lia(lia > 3*num_chn/4) = 1;
    lia = logical(lia);

    loc_event_mat_f(channel,1:length(test(lia))) = test(lia);
    test_peak = peak_event_mat(channel,:);
    test_peak = test_peak(~isnan(test_peak));
    peak_event_mat_f(channel,1:length(test_peak(lia))) = test_peak(lia);
    num_event(channel) = length(test(lia));
end

% Plot Selection Results
figure
subplot(2,1,1)
time = (1/(FS))*(0:length(data)-1/FS);
for chn = 1:num_chn
    clear loc_event_f
    loc_event_f =  loc_event_mat_f(chn,:);
    loc_event_f = loc_event_f(~isnan(loc_event_f));
    plot(loc_event_f,chn*ones(1,length(loc_event_f)),'*');
    hold on 
end
set(gca,'Fontsize',20);
xlim([0 time(end)]);
xlabel('Time (s)');
ylabel('Channel');
title(['Detected Events vs. Channel; Timestamp = ',num2str(data_name)])
set(gca,'YDir','reverse');

% Plot Selection Result of viewing channel
data_abs_view = abs(data(:,view_chn));
subplot(2,1,2)
loc_event_f =  loc_event_mat_f(view_chn,:);
loc_event_f = loc_event_f(~isnan(loc_event_f));
plot(loc_event_f,0.2*ones(1,length(loc_event_f)),'*');
hold on
plot(time,data_abs_view);
set(gca,'Fontsize',20);
xlabel('Time (s)');
ylabel('Amplitude (Pa)');
title(['Event Selection Outcome; Timestamp = ',num2str(data_name),' Channel = ',num2str(view_chn)]);
end
