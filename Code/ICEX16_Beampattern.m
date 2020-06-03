
%% Element Locations
N = 32;
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

%% Plot array Geometry
plot(zeros(length(z),1),z,'o', 'MarkerFaceColor', 'b');
xlabel('x (m)');
ylabel('z (m)');
title('ICEX16 Array Geometry');
set(gca,'Fontsize',30);
grid on

%% Plot Beampattern vs frequency
clear beampattern_mat
%set(0,'DefaultFigureVisible','off')
N = 32;

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/500:pi;
phis = 0;
thetas = pi/2;
weight = 'icex_hanning';
wls = 'none';

f = 2048;
c = 1435;

for i = 1:length(f) 
    disp(i)
    lambda = c/f(i);
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f(i));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %beampattern_mat(i) = getframe(gcf);
end

% v_beampattern = VideoWriter('ICEX16_Array_Beampattern.avi','Uncompressed AVI');
% v_beampattern.FrameRate = 3;
% open(v_beampattern)
% writeVideo(v_beampattern,beampattern_mat)
% close(v_beampattern)

%% Plot Beampattern width
clear beampattern_mat
N = 32;

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/500:pi;
phis = 0;
thetas = pi/2;
weight = 'icex_hanning';
wls = 'none';

f = 850;
c = 1430;

for i = 1:length(f) 
    disp(i)
    lambda = c/f(i);
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f(i));
    [~,loc] = findpeaks(-abs(B));
    angle = -1*(rad2deg(theta_vec)-90);
    half_beamwidth(i) = abs(angle(end)-angle(loc(end)));
    clear loc
end

%% Plot Beampattern vs Elevation
% clear beampattern_mat
% set(0,'DefaultFigureVisible','off')
% N = 32;
% z = zeros(N,1);
% for n = 0:N-1
%     z(n+1) = -(n-(N-1)/2)*d; 
% end

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/100:pi;
phis = 0;
thetas = pi/2;%:pi/36:pi;
weight = 'icex_hanning';
wls = 'none';

f = 900;
c = 1440;
lambda = c/f;

%for i = 1:length(thetas)
   % disp(i)
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %beampattern_mat(i) = getframe(gcf);
%end
% 
% v_beampattern = VideoWriter('ICEX16_Array_Beampattern2.avi','Uncompressed AVI');
% v_beampattern.FrameRate = 3;
% open(v_beampattern)
% writeVideo(v_beampattern,beampattern_mat)
% close(v_beampattern)

