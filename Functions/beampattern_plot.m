
% phi_vec is a vector of angles in phi (in x-y plane).
% theta_vec is a vector angles in theta (measured w/ respect to positive z-axis.
% phis and thetas are the angular angles that the array will be steered to.
% weight is the weighting of the array ('uniform', 'triangluar', etc...)
% p is a matrix that contains the positions of the array elements.
% p = [px1 px2 px3 ... pxN
%      py1 py2 py3 ... pyN 
%      pz1 pz2 pz3 ... pzN]

% Example:
% p = [ zeros(1,20) ; zeros(1,20) ; 9.5:-1:-9.5];
% phi_vec = 0;
% theta_vec = 0:pi/500:pi;
% phis = 0;
% thetas = pi/2;
% weight = 'triangular';
% beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight);
% l/d_ratio = 2;

function [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f)
%d = 0.12;
%lambda = ld_ratio;
N = size(p,2); %number of elements

B = zeros(length(phi_vec),length(theta_vec));
%phi_plot_vec = zeros(length(phi_vec),length(theta_vec));
for ang1 = 1:length(phi_vec)
    phi = phi_vec(ang1);
    for ang2 = 1:length(theta_vec)
        theta = theta_vec(ang2);
        
        k = (-2*pi/lambda)*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)]; % wavenumber
        %phi_plot_vec(ang1,ang2) = -k(3,1)*d;
        ks = (-2*pi/lambda)*[sin(thetas)*cos(phis); sin(thetas)*sin(phis); cos(thetas)]; % wavenumber in steering direction
        
        v = zeros(N,1);
        vs = zeros(N,1);
        for elem = 1:N
            v(elem,1) = exp(-1i*k.'*p(:,elem));
            vs(elem,1) = exp(-1i*ks.'*p(:,elem));
        end


        w = zeros(1,N); 

        if strcmp(weight,'uniform')
            w = (1/N) * ones(1,N);
        end
        
          % icex window
        if strcmp(weight,'icex_hanning')
            w = hanning(42)';
            w(2) = []; w(3) = []; w(4) = []; w(5) = []; w(6) = [];
            w(end-1) = []; w(end-2) = []; w(end-3) = []; w(end-4) = []; w(end-5) = [];
            w = w./sum(w);
        end
        
        % simi xarray window
        if strcmp(weight,'simi_xarray_hanning')
            w = [0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9293 0.7371 0.3694 0.9293  0.7371 0.2248  0.8992 0.7371  0.3162  0.9293 0.7371 0.2594 0.4984 0.6375 0.6341 0.4559 0.6014 0.5083 2.8988e-05 0.4980];  
            w = w./sum(w);
        end
        
        % Hamming window
        if strcmp(weight,'hamming')
            w = hamming(N)';
            w = w./sum(w);
        end
        
        % Hamming window
        if strcmp(weight,'hanning')
            w = hanning(N)';
            w = w./sum(w);
        end
        
        if strcmp(weight,'kaiser')
            w = kaiser(N,2)';
            w = w./sum(w);
        end
        
        if strcmp(weight,'triangular')
            w = triangularPulse(p(3,end), 0, p(3,1), p(3,:));
            w = w/norm(w,1);
        end
        
        if strcmp(weight,'ls')
        w = wls';
        end
        
        wt = w.*vs'; % weight vector transpose
        B(ang1,ang2) = wt*v;

    end
end


% figure
% subplot(1,2,1)
% polarplot(theta_vec,abs(B(1,:))./max(abs(B(1,:))),'linewidth',3);
% thetalim([0 180]);
% title(['Beampattern vs. \theta']);
% set(gca,'Fontsize',30);
% 
% subplot(1,2,2)
% plot(phi_plot_vec(1,:),abs(B(1,:))./max(abs(B(1,:))),'linewidth',3);
% xmin = phi_plot_vec(1,end);
% xmax = phi_plot_vec(1,1);
% grid on
% set(gca,'Fontsize',20);
% title('Beampattern vs. Phi');
% xlabel('\psi_z = (2\picos(\theta)d) / \lambda');
% ylabel('Amplitude');
% xlim([xmin xmax]);


figure
plot(-1*(rad2deg(theta_vec)-90),20*log10(abs(B(1,:))),'linewidth',3);
grid on
set(gca,'Fontsize',30);
title(['Beampattern vs. \theta']);
xlim([0 180]);
%ylim([0 1]);
xlabel('Elevation Angle');
ylabel('Amplitude');
title(['Frequency = ' num2str(f) 'Hz']);

end