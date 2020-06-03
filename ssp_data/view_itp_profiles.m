% view_itp_profiles
% little movie for Rui

% eeshan bhatt
% apr 3, 2020

%% load dataset 
load itp_RUI

%% Plot yearly mean/var

[itp_tt_s,order] = sort(itp_tt);
grid_ssp = grid_ssp(:,order);
itp_lat = itp_lat(order);
itp_lon = itp_lon(order);
itp_date = datestr(itp_tt_s);

figure

subplot(1,2,1);
plot(nanmean(grid_ssp(:,1:1756),2),grid_z,'m-','linewidth',2);
hold on
plot(nanmean(grid_ssp(:,1757:2017),2),grid_z,'c-','linewidth',2);
plot(nanmean(grid_ssp(:,2018:2064),2),grid_z,'b','linewidth',2);
plot(nanmean(grid_ssp(:,2065:2102),2),grid_z,'g-','linewidth',2);
plot(nanmean(grid_ssp(:,2103:2122),2),grid_z,'y-','linewidth',2);
plot(nanmean(grid_ssp(:,2123:2306),2),grid_z,'r-','linewidth',2);
plot(nanmean(grid_ssp(:,2307:end),2),grid_z,'k-','linewidth',2);
plot(nanmean(grid_ssp,2),grid_z,'--','linewidth',2);
grid on
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'ydir','reverse');
set(gca,'Fontsize',20);
title('Mean of Profiles');
legend(['2014 (n = 1756)'],['2015 (n = 261)'],['2016 (n = 47)'],['2017 (n = 38)'],['2018 (n = 20)'],['2019 (n = 184)'],['2020 (n = 26)'],['All (n = 2331)'])

subplot(1,2,2);
plot(sqrt(var(grid_ssp(:,1:1756),0,2,'omitnan')),grid_z,'m-','linewidth',2);
hold on
plot(sqrt(var(grid_ssp(:,1757:2017),0,2,'omitnan')),grid_z,'c-','linewidth',2);
plot(sqrt(var(grid_ssp(:,2018:2064),0,2,'omitnan')),grid_z,'b','linewidth',2);
plot(sqrt(var(grid_ssp(:,2065:2102),0,2,'omitnan')),grid_z,'g-','linewidth',2);
plot(sqrt(var(grid_ssp(:,2103:2122),0,2,'omitnan')),grid_z,'y-','linewidth',2);
plot(sqrt(var(grid_ssp(:,2123:2306),0,2,'omitnan')),grid_z,'r-','linewidth',2);
plot(sqrt(var(grid_ssp(:,2307:end),0,2,'omitnan')),grid_z,'k-','linewidth',2);
plot(sqrt(var(grid_ssp,0,2,'omitnan')),grid_z,'--','linewidth',2);
grid on
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'ydir','reverse');
set(gca,'Fontsize',20);
title('Standard Deviation of Profiles');
legend(['2014 (n = 1756)'],['2015 (n = 261)'],['2016 (n = 47)'],['2017 (n = 38)'],['2018 (n = 20)'],['2019 (n = 184)'],['2020 (n = 26)'],['All (n = 2331)'])


%% Plot peak mean/var
ssp_1peak = zeros(200,length(datetime_1peak));
ssp_2peak = zeros(200,length(datetime_2peak));
ssp_3peak = zeros(200,length(datetime_3peak));

datetimes = datetime(string(itp_date));

for tt = 1:length(datetime_1peak)
    ind = find(datetimes==datetime_1peak(tt));
    ssp_1peak(:,tt) = grid_ssp(:,ind(1));
end

for tt = 1:length(datetime_2peak)
    ind = find(datetimes==datetime_2peak(tt));
    ssp_2peak(:,tt) = grid_ssp(:,ind(1));
end

for tt = 1:length(datetime_3peak)
    ind = find(datetimes==datetime_3peak(tt));
    ssp_3peak(:,tt) = grid_ssp(:,ind(1));
end


figure

subplot(1,2,1);
plot(nanmean(ssp_1peak,2),grid_z,'k-','linewidth',2);
hold on
plot(nanmean(ssp_2peak,2),grid_z,'r-','linewidth',2);
plot(nanmean(ssp_3peak,2),grid_z,'b-','linewidth',2);
grid on
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'ydir','reverse');
set(gca,'Fontsize',20);
title('Mean of Profiles');
legend(['1 peak (n = 1936)'],['2 peaks (n = 205)'],['>3 peaks (n = 70)']);

subplot(1,2,2);
plot(sqrt(var(ssp_1peak,0,2,'omitnan')),grid_z,'k-','linewidth',2);
hold on
plot(sqrt(var(ssp_2peak,0,2,'omitnan')),grid_z,'r-','linewidth',2);
plot(sqrt(var(ssp_3peak,0,2,'omitnan')),grid_z,'b-','linewidth',2);

grid on
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'ydir','reverse');
set(gca,'Fontsize',20);
title('Standard Deviation of Profiles');
legend(['1 peak (n = 1936)'],['2 peaks (n = 205)'],['>3 peaks (n = 70)']);

%% plot peak locationality
figure
hold on
plot(0,0,'k*')
plot(0,0,'r*')
plot(0,0,'b*')
xlabel('Longitude')
ylabel('Latitude')
xlim([-160 -140])
ylim([72.5 75.5])
grid on
set(gca,'FontSize',20)
legend(['1 peak'],['2 peaks'],['>3 peaks'],'AutoUpdate','off');

for jj = 1:length(itp_lat)
    jj
    if ssp_stat(jj).numpeaks==1
        color = 'k';
        plot(itp_lon(jj),itp_lat(jj),[color '*'])
    end
    if ssp_stat(jj).numpeaks==2
        color = 'r';
        plot(itp_lon(jj),itp_lat(jj),[color '*'])
    end
    if ssp_stat(jj).numpeaks>=3
        color = 'b';
        plot(itp_lon(jj),itp_lat(jj),[color '*'])
    end

end

%% Plot peak time
figure
hold on
plot(datetimes(1)-100,0,'k*')
plot(datetimes(1)-100,0,'r*')
plot(datetimes(1)-100,0,'b*')
xlabel('Time')
ylabel('Number of Peaks in SSP')
xlim([datetimes(1) datetimes(end)])
ylim([0 4])
grid on
set(gca,'FontSize',20)
legend(['1 peak'],['2 peaks'],['>3 peaks'],'AutoUpdate','off');

for jj = 1:length(datetimes)
    jj
    if ssp_stat(jj).numpeaks==1
        color = 'k';
        plot(datetimes(jj),ssp_stat(jj).numpeaks,[color '*'])
    end
    if ssp_stat(jj).numpeaks==2
        color = 'r';
        plot(datetimes(jj),ssp_stat(jj).numpeaks,[color '*'])
    end
    if ssp_stat(jj).numpeaks>=3
        color = 'b';
        plot(datetimes(jj),3,[color '*'])
    end
    

end

%% SSP Stat 1 peak
peak_depth_1peak = [];
peak_ssp_1peak = [];
peak_width_1peak = [];
min_depth_1peak = [];
min_ssp_1peak = [];
for ss = 1:length([ssp_stat.numpeaks])
    ss
    if ssp_stat(ss).numpeaks==1
        peak_depth_1peak = [peak_depth_1peak ssp_stat(ss).peakdepths];
        peak_ssp_1peak = [peak_ssp_1peak ssp_stat(ss).peakssps];
        peak_width_1peak = [peak_width_1peak ssp_stat(ss).peak_width];
        min_depth_1peak = [min_depth_1peak min(ssp_stat(ss).mindepth)];
        min_ssp_1peak = [min_ssp_1peak min(ssp_stat(ss).minssp)];
    end
    
end

mean_pd = mean(peak_depth_1peak)
std_pd = sqrt(var(peak_depth_1peak))
mean_pssp = mean(peak_ssp_1peak)
std_pssp = sqrt(var(peak_ssp_1peak))
mean_pw = mean(peak_width_1peak)
std_pw = sqrt(var(peak_width_1peak))
mean_md = mean(min_depth_1peak)
std_md = sqrt(var(min_depth_1peak))
mean_mssp = mean(min_ssp_1peak)
std_mssp = sqrt(var(min_ssp_1peak))


%% oasn

p1mean_nl = reshape(nl3p,56,35);

freqs = [1,3,5,8,10,15,20,25,30,35];
colors = ['k','r','m','c','g','y','b',' ',' ',' '];
figure

subplot(1,4,1)
plot(SSP3p,Depth3p,'linewidth',2);
set(gca,'YDir','reverse');
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)');
grid on
set(gca,'Fontsize',20)
title('3-peak SSP ex.')
hold on
%plot([0 1500],[47.5 47.5],'r--')
%plot([0 1500],[66.7 66.7],'r--')
%plot([0 1500],[78.9 78.9],'r--')
%plot([0 1500],[221.7 221.7],'r--')
xlim([1430 1470]);

subplot (1,4,[2,3,4])
hold on
for i = 1:length(freqs)
    plot(p1mean_nl(:,freqs(i)),nl_depth,[colors(i)],'linewidth',1.5)
end
set(gca,'YDir','reverse');
xlabel('Noise Level (dB)')
yticklabels([]);
grid on
set(gca,'Fontsize',20)
title('Noise Level (35km Monopole Source)')
legend('100Hz','300Hz','500Hz','800Hz','1000Hz','1500Hz','2000Hz','2500Hz','3000Hz','3500Hz','AutoUpdate','off')
hold on
%plot([0 1500],[47.5 47.5],'r--')
%%plot([0 1500],[66.7 66.7],'r--')
%plot([0 1500],[78.9 78.9],'r--')
%plot([0 1500],[221.7 221.7],'r--')
xlim([20 70])

%% oasn discrete range

mssp05km = dlmread('meanssp_20km_new.mtv',' ');
mssp05km = mssp05km(1:1960,2);
mssp05km = reshape(mssp05km,56,35);

nl_depth = [5:10:305 325 340 360 380 400 420 440 460 480 500:20:800];
% simi94 = dlmread('simi94_new.mtv',' ');
% Depth = simi94(1:79,2);
% SSP = simi94(1:79,1);
% siminl = simi94(89:2048,2);
% siminl = reshape(siminl,56,35);

freqs = [1,3,5,8,10,15,20,25,30,35];
colors = ['k','r','m','c','g','y','b',' ',' ',' '];
figure

subplot(1,4,1)
plot(SSP,Depth,'linewidth',2);
set(gca,'YDir','reverse');
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)');
grid on
set(gca,'Fontsize',20)
title('1-Peak Mean SSP')
hold on
plot([0 1500],[65 65],'r--','linewidth',1.5)
plot([0 1500],[213.7 213.7],'r--','linewidth',1.5)
xlim([1430 1470]);
ylim([0 800]);

subplot (1,4,[2,3,4])
hold on
for i = 1:length(freqs)
    plot(mssp05km(:,freqs(i)),nl_depth,[colors(i)],'linewidth',1.5)
end
set(gca,'YDir','reverse');
xlabel('Noise Level (dB)')
yticklabels([]);
grid on
set(gca,'Fontsize',20)
title('Noise Level (20Km Discrete Source)')
legend('100Hz','300Hz','500Hz','800Hz','1000Hz','1500Hz','2000Hz','2500Hz','3000Hz','3500Hz','AutoUpdate','off')
hold on
plot([0 1500],[65 65],'r--','linewidth',1.5)
plot([0 1500],[213.7 213.7],'r--','linewidth',1.5)
xlim([0 80])
ylim([0 800]);

%% oasn discrete range freq

mssp0km = dlmread('p2ssp2_0km_new.mtv',' ');
mssp0km = mssp0km(1:1960,2);
mssp0km = reshape(mssp0km,56,35);

% mssp05km = dlmread('p2ssp2_0.5km_new.mtv',' ');
% mssp05km = mssp05km(1:1960,2);
% mssp05km = reshape(mssp05km,56,35);

mssp1km = dlmread('p2ssp2_1km_new.mtv',' ');
mssp1km = mssp1km(1:1960,2);
mssp1km = reshape(mssp1km,56,35);

% mssp15km = dlmread('meanssp_1.5km_new.mtv',' ');
% mssp15km = mssp15km(1:1960,2);
% mssp15km = reshape(mssp15km,56,35);

mssp2km = dlmread('p2ssp2_2km_new.mtv',' ');
mssp2km = mssp2km(1:1960,2);
mssp2km = reshape(mssp2km,56,35);

% mssp25km = dlmread('meanssp_2.5km_new.mtv',' ');
% mssp25km = mssp25km(1:1960,2);
% mssp25km = reshape(mssp25km,56,35);

mssp3km = dlmread('p2ssp2_3km_new.mtv',' ');
mssp3km = mssp3km(1:1960,2);
mssp3km = reshape(mssp3km,56,35);

% mssp35km = dlmread('meanssp_3.5km_new.mtv',' ');
% mssp35km = mssp35km(1:1960,2);
% mssp35km = reshape(mssp35km,56,35);

mssp4km = dlmread('p2ssp2_4km_new.mtv',' ');
mssp4km = mssp4km(1:1960,2);
mssp4km = reshape(mssp4km,56,35);

% mssp45km = dlmread('meanssp_4.5km_new.mtv',' ');
% mssp45km = mssp45km(1:1960,2);
% mssp45km = reshape(mssp45km,56,35);

mssp5km = dlmread('p2ssp2_5km_new.mtv',' ');
mssp5km = mssp5km(1:1960,2);
mssp5km = reshape(mssp5km,56,35);

mssp6km = dlmread('p2ssp2_6km_new.mtv',' ');
mssp6km = mssp6km(1:1960,2);
mssp6km = reshape(mssp6km,56,35);

mssp7km = dlmread('p2ssp2_7km_new.mtv',' ');
mssp7km = mssp7km(1:1960,2);
mssp7km = reshape(mssp7km,56,35);

mssp8km = dlmread('p2ssp2_8km_new.mtv',' ');
mssp8km = mssp8km(1:1960,2);
mssp8km = reshape(mssp8km,56,35);

mssp9km = dlmread('p2ssp2_9km_new.mtv',' ');
mssp9km = mssp9km(1:1960,2);
mssp9km = reshape(mssp9km,56,35);

mssp10km = dlmread('p2ssp2_10km_new.mtv',' ');
mssp10km = mssp10km(1:1960,2);
mssp10km = reshape(mssp10km,56,35);

mssp105km = dlmread('p2ssp2_15km_new.mtv',' ');
mssp105km = mssp105km(1:1960,2);
mssp105km = reshape(mssp105km,56,35);

mssp20km = dlmread('p2ssp2_20km_new.mtv',' ');
mssp20km = mssp20km(1:1960,2);
mssp20km = reshape(mssp20km,56,35);

mssp205km = dlmread('p2ssp2_25km_new.mtv',' ');
mssp205km = mssp205km(1:1960,2);
mssp205km = reshape(mssp205km,56,35);

mssp30km = dlmread('p2ssp2_30km_new.mtv',' ');
mssp30km = mssp30km(1:1960,2);
mssp30km = reshape(mssp30km,56,35);

f = 30;
dist = [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 9 10 15 20 25 30];
colors = ['k','r','m','c','g','y','b',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '];
figure

subplot(1,4,1)
plot(SSP,Depth,'linewidth',2);
set(gca,'YDir','reverse');
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)');
grid on
set(gca,'Fontsize',20)
title('2-Peak SSP')
hold on
plot([0 1500],[54.9 54.9],'r--')
plot([0 1500],[105.9 105.9],'r--')
plot([0 1500],[256 256],'r--')
xlim([1430 1470]);

subplot (1,4,[2,3,4])
hold on

%plot(mssp0km(:,f),nl_depth,[colors(1)],'linewidth',1.5)
%plot(mssp05km(:,f),nl_depth,[colors(1)],'linewidth',1.5)
plot(mssp1km(:,f),nl_depth,[colors(1)],'linewidth',1.5)
%plot(mssp15km(:,f),nl_depth,[colors(4)],'linewidth',1.5)
%plot(mssp2km(:,f),nl_depth,[colors(3)],'linewidth',1.5)
%plot(mssp25km(:,f),nl_depth,[colors(6)],'linewidth',1.5)
plot(mssp3km(:,f),nl_depth,[colors(2)],'linewidth',1.5)
%plot(mssp35km(:,f),nl_depth,[colors(2)],'linewidth',1.5)
%plot(mssp4km(:,f),nl_depth,[colors(9)],'linewidth',1.5)
%plot(mssp45km(:,f),nl_depth,[colors(10)],'linewidth',1.5)
%plot(mssp5km(:,f),nl_depth,[colors(4)],'linewidth',1.5)
plot(mssp6km(:,f),nl_depth,[colors(3)],'linewidth',1.5)
%plot(mssp7km(:,f),nl_depth,[colors(13)],'linewidth',1.5)
%plot(mssp8km(:,f),nl_depth,[colors(14)],'linewidth',1.5)
%plot(mssp9km(:,f),nl_depth,[colors(15)],'linewidth',1.5)
plot(mssp10km(:,f),nl_depth,[colors(4)],'linewidth',1.5)
plot(mssp105km(:,f),nl_depth,[colors(5)],'linewidth',1.5)
plot(mssp20km(:,f),nl_depth,[colors(6)],'linewidth',1.5)
plot(mssp205km(:,f),nl_depth,[colors(7)],'linewidth',1.5)
%plot(mssp30km(:,f),nl_depth,[colors(20)],'linewidth',1.5)

set(gca,'YDir','reverse');
xlabel('Noise Level (dB)')
yticklabels([]);
grid on
set(gca,'Fontsize',20)
title('Noise Level (Frequency = 3000 Hz)')
legend('1Km','3Km','6Km','10Km','15Km','20Km','25Km','AutoUpdate','off')
hold on
plot([0 1500],[54.9 54.9],'r--')
plot([0 1500],[105.9 105.9],'r--')
plot([0 1500],[256 256],'r--')
xlim([0 100])

%% se sim
freq = 100:100:3500;
range = [10 1000:1000:50000];

sl = dlmread('source_level_new.mtv',' ');
sl = sl(1:1785,2);
sl = reshape(sl,51,35).';

unl = dlmread('uniform_nl_new.mtv',' ');
unl = unl(1:1785,2);
unl = reshape(unl,51,35).';


dnl = dlmread('discrete_nl_new.mtv',' ');
dnl = dnl(1:1785,2);
dnl = reshape(dnl,51,35).';

%dnl(:,ceil(end/2):end) = dnl(:,1:ceil(end/2));
%dnl(:,1:floor(end/2)) = fliplr(dnl(:,ceil(end/2)+1:end));

usnr = sl-unl(:,1);
dsnr = sl-dnl;

figure
h1 = pcolor(range,freq,usnr);
set(h1,'Edgecolor','None');
xlabel('Range (m)');
ylabel('Frequency (Hz)')
set(gca,'Fontsize',20);
caxis([0 160])
colorbar
colormap jet

figure
h2 = pcolor(range,freq,dsnr);
set(h2,'Edgecolor','None');
xlabel('Range (m)');
ylabel('Frequency (Hz)')
set(gca,'Fontsize',20);
caxis([0 160])
colorbar
colormap jet


%% figure
figure(1); clf

% mean
plot(mean(grid_ssp,2,'omitnan'),grid_z,'r--');
set(gca,'ydir','reverse')
grid on
xlabel('ssp [m/s]')
ylabel('depth [m]')

% animated line for profile
an_profile = animatedline('linewidth',4,'color','b');

for k = 1:length(itp_ssp)
    clearpoints(an_profile);
    addpoints(an_profile,itp_ssp{k},itp_depth{k});
    title(datestr(itp_tt(k)));
    drawnow
end

%% Find Peaks

[itp_tt_s,order] = sort(itp_tt);
grid_ssp = grid_ssp(:,order);

datetime_1peak = [];
datetime_2peak = [];
datetime_3peak = [];

for k = 1:length(itp_ssp)
    k
    temp = grid_ssp(15:70,k);
    temp1 = grid_ssp(:,k);
    
    [pks,locs,w,p] = findpeaks(temp,grid_z(15:70),'MinPeakProminence',0.25);
    ssp_stat(k).numpeaks = length(pks);
    if length(pks) == 1
        datetime_1peak = [datetime_1peak datetimes(k)];
    end
    
    if length(pks) == 2
        datetime_2peak = [datetime_2peak datetimes(k)];
    end
    
    if length(pks) >= 3
        datetime_3peak = [datetime_3peak datetimes(k)];
    end
        
    ssp_stat(k).peakdepths = locs;
    ssp_stat(k).peakssps = pks;
    ssp_stat(k).peak_width = 2*w;
    ssp_stat(k).peak_prominence = p;
    
    tf = islocalmin(temp1,'MinProminence',1.25);
    ssp_stat(k).mindepth = grid_z(tf);
    ssp_stat(k).minssp = temp1(tf);
    
end

%% Peak stats

[itp_tt_s,order] = sort(itp_tt);
grid_ssp = grid_ssp(:,order);

for kk = 1:length(ssp_stat)
    
    if isempty(ssp_stat(kk).peakdepths)
        continue
        
    elseif length(ssp_stat(kk).peakdepths) == 1
        %1peak_cases(kk) = ssp_stat(kk);
        onepeak_cases(kk).time = itp_tt(kk);
        onepeak_cases(kk).mppeak_depth = ssp_stat(kk).peakdepths;
        onepeak_cases(kk).mppeak_ssp = ssp_stat(kk).peakssps;
        onepeak_cases(kk).mppeak_halfwidth = ssp_stat(kk).peak_halfwidth;
        onepeak_cases(kk).mppeak_peakprom = ssp_stat(kk).peak_prominence;
        onepeak_cases(kk).mmpeak_ssp = ssp_stat(kk).minssp;
        onepeak_cases(kk).mmpeak_depth = ssp_stat(kk).mindepth;
        
    elseif length(ssp_stat(kk).peakdepths) == 1
        %2peak_cases(kk) = ssp_stat(kk);
        twopeak_cases(kk).time = itp_tt(kk);
        [~,num] = max(ssp_stat(kk).peak_prominence);
        twopeak_cases(kk).mppeak_depth = ssp_stat(kk).peakdepths(num);
        twopeak_cases(kk).mppeak_ssp = ssp_stat(kk).peakssps(num);
        twopeak_cases(kk).mppeak_halfwidth = ssp_stat(kk).peak_halfwidth(num);
        twopeak_cases(kk).mppeak_peakprom = ssp_stat(kk).peak_prominence(num);
        [val,num2] = min(ssp_stat(kk).minssp);
        twopeak_cases(kk).mmpeak_ssp = val;
        twopeak_cases(kk).mmpeak_depth = ssp_stat(kk).mindepth(num2);
        
    else
       %3peak_cases(kk) = ssp_stat(kk);
       threepeak_cases(kk).time = itp_tt(kk);
       [~,num] = max(ssp_stat(kk).peak_prominence);
       threepeak_cases(kk).mppeak_depth = ssp_stat(kk).peakdepths(num);
       threepeak_cases(kk).mppeak_ssp = ssp_stat(kk).peakssps(num);
       threepeak_cases(kk).mppeak_halfwidth = ssp_stat(kk).peak_halfwidth(num);
       threepeak_cases(kk).mppeak_peakprom = ssp_stat(kk).peak_prominence(num);
       [val,num2] = min(ssp_stat(kk).minssp);
       threepeak_cases(kk).mmpeak_ssp = val;
       threepeak_cases(kk).mmpeak_depth = ssp_stat(kk).mindepth(num2);
           
    end 
      
        
     
    
    
    
end
