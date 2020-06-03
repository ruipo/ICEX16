

gaplength_mat = [];
gaplength_g = [];
clusterlength_mat = [];
clusternum_mat = [];
clusterrate_mat = [];
clusterlength_g = [];

%%
a = 6;
TimeBetweenEvents = TimeBetweenEvents6(~isnan(TimeBetweenEvents6));
clear timelog

timelog(1) = 0;
for i = 1:length(TimeBetweenEvents)
    timelog(i+1) = timelog(i) + TimeBetweenEvents(i);
end
%timelog(length(timelog)+1) = 480*60;

figure
h = histogram(timelog,480,'BinLimits',[0 28800]); %SIMI94
%h = histogram(timelog,339,'BinLimits',[0 20340]); %ICEX16
xlabel('Time (Hours)')
ylabel('Number of Events');
title('SIMI94 80-160 Hz');
set(gca,'Fontsize',30);
eventcluster = h.Values;

% determine gap lengths
k = find(eventcluster);
for i = 1:length(k)-1
    gaplength(i) = k(i+1)-k(i) - 1;
end
gaplength = gaplength(gaplength > 0);

if k(1) ~= 1;
    gaplength = [k(1)-1 gaplength]; % length of gaps between clusters (in min)
end

mean_gaplength(a) = mean(gaplength);
gaplength_mat = [gaplength_mat gaplength];
gaplength_g = [gaplength_g a*ones(size(gaplength))];


% determine cluster lengths and number of eventer in clusters
c = find(~eventcluster);
for i = 1:length(c)-1
    clusterlength(i) = c(i+1)-c(i) - 1;
    clusternum(i) = sum(eventcluster(c(i):c(i+1))); 
end
clusterlength = clusterlength(clusterlength > 0);
clusternum = clusternum(clusternum > 0);

if c(1) ~= 1;
    clusterlength = [c(1)-1 clusterlength]; % lengths of event clusters (in min)
    clusternum = [eventcluster(1) clusternum]; % number of events in each cluster
end

clusterrate = clusternum./clusterlength; % events per min in each cluster

mean_clusterlength(a) = mean(clusterlength);
mean_clusternum(a) = mean(clusternum);
mean_clusterrate(a) = mean(clusterrate);
clusterlength_mat = [clusterlength_mat clusterlength];
clusternum_mat = [clusternum_mat clusternum];
clusterrate_mat = [clusterrate_mat clusterlength];
clusterlength_g = [clusterlength_g a*ones(size(clusterlength))];


%%

figure
h1 = boxplot(gaplength_mat,gaplength_g);
ylabel('Gaps between Event Clusters (min)');
set(gca,'XTickLabel',{'ICEX16','SIMI94','ICEX16','SIMI94','ICEX16','SIMI94'})
set(gca,'Fontsize',30);
set(h1,{'linew'},{2})
set(gca, 'YScale', 'log')
ylim([0.5 100]);
grid on

figure
h2 = boxplot(clusterlength_mat,clusterlength_g);
ylabel('Length of Event Clusters (min)');
set(gca,'XTickLabel',{'ICEX16','SIMI94','ICEX16','SIMI94','ICEX16','SIMI94'})
set(gca,'Fontsize',30);
set(h2,{'linew'},{2})
set(gca, 'YScale', 'log')
ylim([0.5 100]);
grid on

figure
h3 = boxplot(clusternum_mat,clusterlength_g);
ylabel('Number of Events in Each Event Cluster');
set(gca,'XTickLabel',{'ICEX16','SIMI94','ICEX16','SIMI94','ICEX16','SIMI94'})
set(gca,'Fontsize',30);
set(h3,{'linew'},{2})
set(gca, 'YScale', 'log')
ylim([0.5 5000]);
grid on

figure
h4 = boxplot(clusterrate_mat,clusterlength_g);
ylabel('Average Event Rate in Each Cluster (events/min)');
set(gca,'XTickLabel',{'ICEX16','SIMI94','ICEX16','SIMI94','ICEX16','SIMI94'})
set(gca,'Fontsize',30);
set(h4,{'linew'},{2})
set(gca, 'YScale', 'log')
ylim([0.5 100]);
grid on

%%


SIMI_mean_clusterlength = [mean_clusterlength(1) mean_clusterlength(3) mean_clusterlength(5)];
SIMI_mean_clusternum = [mean_clusternum(1) mean_clusternum(3) mean_clusternum(5)];
SIMI_mean_clusterrate = [mean_clusterrate(1) mean_clusterrate(3) mean_clusterrate(5)];
SIMI_mean_gaplength = [mean_gaplength(1) mean_gaplength(3) mean_gaplength(5)];
ICEX_mean_gaplength = [mean_gaplength(2) mean_gaplength(4) mean_gaplength(6) mean_gaplength(7) mean_gaplength(8)];
ICEX_mean_clusterlength = [mean_clusterlength(2) mean_clusterlength(4) mean_clusterlength(6) mean_clusterlength(7) mean_clusterlength(8)];
ICEX_mean_clusternum = [mean_clusternum(2) mean_clusternum(4) mean_clusternum(6) mean_clusternum(7) mean_clusternum(8)];
ICEX_mean_clusterrate = [mean_clusterrate(2) mean_clusterrate(4) mean_clusterrate(6) mean_clusterrate(7) mean_clusterrate(8)];
%%
figure
plot(SIMI_mean_gaplength)
hold on
plot(ICEX_mean_gaplength)
ylabel('Mean Gap Length between Event Clusters (min)');
%set(gca,'XTickLabel',{'40-80 Hz', '80-160 Hz', '160-320 Hz', '320-640 Hz', '640-1280 Hz'})
set(gca,'Fontsize',30);
grid on

figure
plot(SIMI_mean_clusterlength)
hold on
plot(ICEX_mean_clusterlength)
ylabel('Mean Length of Event Clusters (min)');
%set(gca,'XTickLabel',{'40-80 Hz', '80-160 Hz', '160-320 Hz', '320-640 Hz', '640-1280 Hz'})
set(gca,'Fontsize',30);
grid on

figure
plot(SIMI_mean_clusternum)
hold on
plot(ICEX_mean_clusternum)
ylabel('Average Number of Events in Each Event Cluster');
%set(gca,'XTickLabel',{'40-80 Hz', '80-160 Hz', '160-320 Hz', '320-640 Hz', '640-1280 Hz'})
set(gca,'Fontsize',30);
grid on

figure
plot(SIMI_mean_clusterrate)
hold on
plot(ICEX_mean_clusterrate)
ylabel('Average Event Rate in Each Cluster (events/min)');
%set(gca,'XTickLabel',{'40-80 Hz', '80-160 Hz', '160-320 Hz', '320-640 Hz', '640-1280 Hz'})
set(gca,'Fontsize',30);
grid on


