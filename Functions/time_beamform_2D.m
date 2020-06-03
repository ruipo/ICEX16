% INPUTS
%
% data		time series (time x element)
% p		element positions matrix,  [x_e1,y_e1;x_e2...]
% FS		sample frequency (Hz)
% azlist		desired look angles in azmuth in radians
% c		sound speed
% overlap    percent overlap [0 1]
% win_len    window length for data segments
%
% OUTPUTS
%
% time_beamform	   time domain beamformed result (az x time)

function[time_beamform,t] = time_beamform_2D(data,p,FS,azlist,c,win_len,overlap)

N = size(data,2);
window_start = round(win_len-win_len*overlap);
num_window = round(size(data,1)/window_start)-3;

D = zeros(size(p));
t = zeros(num_window-1,1);
dr = zeros(2,N);
dvec = zeros(N,1);
dt = zeros(N,1);
data_seg = zeros(N,win_len);
time_beamform = zeros(length(azlist),num_window-1);
num_window

for i = 1:N
    D(:,i) = p(:,1) - p(:,i);
end

for l = 2:num_window
    l
    t(l-1) = ((l)*window_start-window_start+1)/FS;
    
    for az = 1:length(azlist);
        x = round(cos(azlist(az)),10);
        y = round(sin(azlist(az)),10);

        k = -[x;y];

        for elem = 1:N
            dr(:,elem) = (dot(D(:,elem),k)/(norm(k)^2))*k;
            dvec(elem) = norm(dr(:,elem));
            
            if azlist(az) == 0;
                if dr(1,elem) < 0
                    dvec(elem) = -dvec(elem);
                end
            end

            if 0 < azlist(az) && azlist(az) < pi
                if dr(2,elem) < 0
                    dvec(elem) = -dvec(elem);
                end
            end
            
            if azlist(az) == pi;
                if dr(1,elem) > 0
                    dvec(elem) = -dvec(elem);
                end
            end
            
            if pi < azlist(az)
                if dr(2,elem) > 0
                    dvec(elem) = -dvec(elem);
                end
            end

            dt(elem) = dvec(elem)./c;

            %[data_seg(elem,:),~] = envelope(data(l*window_start-window_start+1+round(dt(elem)*FS):l*window_start-window_start+round(dt(elem)*FS)+win_len,elem));
            data_seg(elem,:) = data(l*window_start-window_start+1+round(dt(elem)*FS):l*window_start-window_start+round(dt(elem)*FS)+win_len,elem);

            ht = hilbert(mean(data_seg));
            upper = abs(ht);
            
            time_beamform(az,l-1) = mean(upper);

        end
              
    end
     
end