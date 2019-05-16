
f=60;
fs=8192;
N=762520;
A=1;
fase=pi/2;

discreto=1:N;
t=discreto/fs;
g1=A*sin(2*pi*f*t);
g2=A*sin(2*pi*f*t+fase);

plot(t,g1);
hold on;
plot(t,g2);
xlim([17.45 17.5]);