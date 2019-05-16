
clear all;
close all;
load senal_ECG1.mat
format long

notch=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'DesignMethod','butter','SampleRate',8192);


A=1;
SNR=20;
SigmaRuido= sqrt (A^2 / (2*10^(SNR/10)));
fs=8192;
frecuencia1=60;
frecuencia2=70;
N=length(ECG_signal1);
fase=pi/2;
Nfiltro=100;

discreto1=1:N/2;
t1=discreto1'/fs;
g1=A*sin(2*pi*frecuencia1*t1+fase);
discreto2=N/2+1:N;
t2=discreto2'/fs;
g2=A*sin(2*pi*frecuencia2*t2+fase);
g_2=[g1;g2];
t=[t1;t2];



discreto=1:N;
t=discreto'/fs;
g=A*sin(2*pi*frecuencia1*t+fase);

d1=ECG_signal1'+g+SigmaRuido*randn;%señal con ruido y con tono de 60
d2=ECG_signal1'+g_2+SigmaRuido*randn;%señal con ruido y con tono de 60/70
d3=ECG_signal1'+g;%señal sin ruido y con tono de 60
d4=ECG_signal1'+g_2;%señal sin ruido y con tono de 60/70

SenalLimpia1=filtfilt(notch,d1);
SenalLimpia2=filtfilt(notch,d2);
SenalLimpia3=filtfilt(notch,d3);
SenalLimpia4=filtfilt(notch,d4);

figure
subplot(2,1,1);
plot(SenalLimpia1);
hold on
plot(ECG_signal1);
xlim([1 6000]);
ylim([-0.5,1.25]);
grid on;
title('Señales de ECG SNR=20 dB f_T=60 Hz');
subplot(2,1,2);
plot(SenalLimpia3);
hold on
plot(ECG_signal1);
title('Señales de ECG  sin ruido f_T=60 Hz');
xlim([1 6000]);
ylim([-0.5,1.25]);
grid on;
% print('ejercicio4c.png','-dpng');

figure
subplot(2,1,1);
plot(SenalLimpia2);
hold on
plot(ECG_signal1);
grid on;
xlim([3.75 3.95]*10^5);
ylim([-1.2 2]);
title('Señales de ECG SNR=20 dB f_T=60/70 Hz');
subplot(2,1,2);
plot(SenalLimpia4);
hold on
plot(ECG_signal1);
grid on;
xlim([3.75 3.95]*10^5);
ylim([-1.2,2]);
title('Señales de ECG sin ruido f_T=60/70 Hz');
% print('ejercicio4c_70.png','-dpng');