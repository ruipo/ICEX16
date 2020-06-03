

function[beamform_output] = beamform_covmat(cov_mat,p,elev,az,c,f,weighting)

N = size(p,1);
beam_elev = (90-elev).*(pi/180);
beam_az = az.*(pi/180);

% Start beamforming

% calculate k
k = 2*pi*f./c;

% linear window
if strcmp(weighting,'uniform')
    win = ones(1,N);
    win = win./norm(win,2);
end

% icex window
if strcmp(weighting,'icex_hanning')
    win = hanning(42)';
    win(2) = []; win(3) = []; win(4) = []; win(5) = []; win(6) = [];
    win(end-1) = []; win(end-2) = []; win(end-3) = []; win(end-4) = []; win(end-5) = [];
    win = win./norm(win,2);
end

% simi xarray window
if strcmp(weighting,'simi_xarray_hanning')
    win = [0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9293 0.7371 0.3694 0.9293  0.7371 0.2248  0.8992 0.7371  0.3162  0.9293 0.7371 0.2594 0.4984 0.6375 0.6341 0.4559 0.6014 0.5083 2.8988e-05 0.4980];  
    win = win./norm(win,2);
end

% Hamming window
 if strcmp(weighting,'hanning')
    win = hanning(N)';
    win = win./norm(win,2);
end

% build steering vectors
for j = 1:length(beam_az)
    for mm = 1:length(beam_elev)
        mm
        % form steering vector
        steer = exp(1i * k *(sin(beam_elev(mm))*cos(beam_az(j))*p(:,1)'+sin(beam_elev(mm))*sin(beam_az(j))*p(:,2)'+cos(beam_elev(mm))*p(:,3)'));

        % apply weighting
        steer = (steer.*(ones(size(k,1),1)*win)).';

        % beamform
        for l = 1:size(cov_mat,3)
            b_elem = steer'*cov_mat(:,:,l)*steer;
            beamform_output(mm,l,j) = abs(b_elem).^2;
        end
    end
end

end