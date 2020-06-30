h = findobj(gca,'Type','line');
x3=get(h,'Xdata');
y3=get(h,'Ydata');

%%

figure

subplot(3,4,[1 5 9])
plot(SSP,Depth,'-','linewidth',2)
ylim([0 3000])
set(gca,'YDir','reverse')
grid on
set(gca,'fontsize',20)
hold on
plot(1445,38 ,'ro','MarkerFaceColor','red');
plot(1445,138 ,'ro','MarkerFaceColor','red');
plot(1445,238 ,'ro','MarkerFaceColor','red');
title('ICEX16 SSP')
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')

subplot(3,4,[2 3 4])
set(gca,'YDir','reverse');
hold on
for i = 1:length(x1)
    plot(x1{i},y1{i},'g','linewidth',1.5)
end
grid on
ylim([0 3000])
xlim([0 50000])
set(gca,'fontsize',20)
title('ICEX16 38m [-8 -15] Deg. Rays')

% xlabel('Range (km)')
xticks([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
xticklabels({'0','5','10','15','20','25','30','35','40','45','50'})

subplot(3,4,[6 7 8])
set(gca,'YDir','reverse');
hold on
for i = 1:length(x2)
    plot(x2{i},y2{i},'g','linewidth',1.5)
end
grid on
ylim([0 3000])
xlim([0 50000])
set(gca,'fontsize',20)
title('ICEX16 138m [-8 -15] Deg. Rays')

% xlabel('Range (km)')
xticks([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
xticklabels({'0','5','10','15','20','25','30','35','40','45','50'})

subplot(3,4,[10 11 12])
set(gca,'YDir','reverse');
hold on
for i = 1:length(x3)
    plot(x3{i},y3{i},'g','linewidth',1.5)
end
grid on
ylim([0 3000])
xlim([0 50000])
set(gca,'fontsize',20)
title('ICEX16 238m [-8 -15] Deg. Rays')

xlabel('Range (km)')
xticks([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
xticklabels({'0','5','10','15','20','25','30','35','40','45','50'})