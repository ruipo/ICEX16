[A,B,C,D] = butter(10,[500 560]/750);

sos = ss2sos(A,B,C,D);

[z,p,k] = butter(50,[20 900]/6000,'bandpass');
sos = zp2sos(z,p,k);
hold on
fvt1 = fvtool(sos,'Fs',12000);

[b,a] = butter(6,[20 1500]/6000,'bandpass');
freqz(b,a);


[z,p,k] = butter(100,[20 900]/6000,'bandpass');
sos = zp2sos(z,p,k);
fvt1 = fvtool(sos,'Fs',12000);
aco_fil = sosfilt(sos,aco_norm(:,16));
aoc_fil = aco_fil';

%%

theta = -pi:0.01:pi;
%N = 16;
N = 22;

l = 1435/f;
%d = 7.5;
d = 0.75;
s = deg2rad(40);

f = 500;
B = 1/N*(sin(N/2*2*pi/l*(cos(theta)-cos(s))*d)./sin(0.5*2*pi/l*(cos(theta)-cos(s))*d));
%BC = 0.5*sin(pi/(2*N))*((sin(N*pi/2*d*(cos(theta)-1/N))./sin(pi/2*d*(cos(theta)-1/N))) + (sin(N*pi/2*d*(cos(theta)+1/N))./sin(pi/2*d*(cos(theta)+1/N))));

polarplot(theta,B)
figure
plot(theta,B)
%polarplot(theta,BC)
%hold on
