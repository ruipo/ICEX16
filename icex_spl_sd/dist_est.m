
%% ICEX16 38m
% get beamform_output_3_50_norm from ICEX_oasn
elev = -90:1:90;
dist = 3:0.5:50;

beamform_output_3_50_norm = zeros(size(beamform_output_3_50));
beamform_mat_norm = zeros(size(beamform_mat(1:end-48,:)));

for i=1:95
    beamform_output_3_50_norm(:,i) = beamform_output_3_50(:,i)./max(beamform_output_3_50(:,i));
end

for i=1:1070-48
    beamform_mat_norm(i,:) = beamform_mat(i,:)./max(beamform_mat(i,:));
end

dist_est = [];
for k1 = 1:size(beamform_mat_norm,1)
    k1
    D = [];
    for k = 1:size(beamform_output_3_50_norm,2)
        D(k) = sqrt(sum((beamform_mat_norm(k1,:) - beamform_output_3_50_norm(:,k).').^2));  
    end
    [~,ind] = min(D);
    dist_est(k1) = dist(ind);
end

figure
subplot(2,1,1)
h = pcolor(time_vec(1:end-48),elev,10*log10(beamform_mat_norm.'));
set(h,'Edgecolor','none')
xticklabels([])
ylabel('Elevation (Degrees)')
colormap jet
caxis([-20 10])
ylim([-25 25])
colorbar
set(gca,'fontsize',20)
title('ICEX16 38m Beamforming output')

subplot(2,1,2)
plot(datetimevec(1:end-48),dist_est,'bo')
set(gca,'fontsize',20)
xlim([datetimevec(1) datetimevec(end-48)])
ylim([3 50])
ylabel('Distance (km)')
xlabel('Time (UTC)')
title('Source Range Estimation using Discrete Source Model')
grid on

%% ICEX16 138m

beamform_output_3_50_norm = zeros(size(beamform_output_3_50));
beamform_mat_norm = zeros(size(beamform_mat(139:160,:)));

for i=1:95
    beamform_output_3_50_norm(:,i) = beamform_output_3_50(:,i)./max(beamform_output_3_50(:,i));
end

for i=1:size(beamform_mat_norm,1)
    beamform_mat_norm(i,:) = beamform_mat(138+i,:)./max(beamform_mat(138+i,:));
end

dist_est = [];
for k1 = 1:size(beamform_mat_norm,1)
    k1
    D = [];
    for k = 1:size(beamform_output_3_50_norm,2)
        D(k) = sqrt(sum((beamform_mat_norm(k1,:) - beamform_output_3_50_norm(:,k).').^2));  
    end
    [~,ind] = min(D);
    dist_est(k1) = dist(ind);
end

figure
subplot(2,1,1)
h = pcolor([0:30:size(beamform_mat_norm,1)*30-30],elev,10*log10(beamform_mat_norm.'));
set(h,'Edgecolor','none')
xticklabels([])
ylabel('Elevation (Degrees)')
colormap jet
caxis([-20 10])
ylim([-25 25])
colorbar
set(gca,'fontsize',20)
title('ICEX16 138m Beamforming output')

subplot(2,1,2)
plot([0:30:size(beamform_mat_norm,1)*30-30],dist_est,'bo')
set(gca,'fontsize',20)
xlim([0 size(beamform_mat_norm,1)*30-30])
ylim([3 50])
ylabel('Distance (km)')
xlabel('Time (Seconds)')
title('Source Range Estimation using Discrete Source Model')
grid on

%% ICEX16 238m

beamform_output_3_50_norm = zeros(size(beamform_output_3_50));
beamform_mat_norm = zeros(size(beamform_mat(300:341,:)));

for i=1:95
    beamform_output_3_50_norm(:,i) = beamform_output_3_50(:,i)./max(beamform_output_3_50(:,i));
end

for i=1:size(beamform_mat_norm,1)
    beamform_mat_norm(i,:) = beamform_mat(299+i,:)./max(beamform_mat(299+i,:));
end

dist_est = [];
for k1 = 1:size(beamform_mat_norm,1)
    k1
    D = [];
    for k = 1:size(beamform_output_3_50_norm,2)
        D(k) = sqrt(sum((beamform_mat_norm(k1,:) - beamform_output_3_50_norm(:,k).').^2));  
    end
    [~,ind] = min(D);
    dist_est(k1) = dist(ind);
end

figure
subplot(2,1,1)
h = pcolor([0:30:size(beamform_mat_norm,1)*30-30],elev,10*log10(beamform_mat_norm.'));
set(h,'Edgecolor','none')
xticklabels([])
ylabel('Elevation (Degrees)')
colormap jet
caxis([-20 10])
ylim([-25 25])
colorbar
set(gca,'fontsize',20)
title('ICEX16 238m Beamforming output')

subplot(2,1,2)
plot([0:30:size(beamform_mat_norm,1)*30-30],dist_est,'bo')
set(gca,'fontsize',20)
xlim([0 size(beamform_mat_norm,1)*30-30])
ylim([3 50])
ylabel('Distance (km)')
xlabel('Time (Seconds)')
title('Source Range Estimation using Discrete Source Model')
grid on