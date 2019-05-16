clear all;
close all;
load senal_ECG1.mat
format long

A=1;
B=0.1*A;
SNR=20;
SigmaRuido= sqrt (A^2 / (2*10^(SNR/10)));
fs=8192;
frecuencia1=60;
frecuencia2=70;
N=length(ECG_signal1);
fase=pi/2;
Nfiltro=100;
mu=0.001;
SenalEstimada=zeros(size(ECG_signal1))';

discreto1=1:N/2;
t1=discreto1'/fs;
g1=A*sin(2*pi*frecuencia1*t1+fase);
discreto2=N/2+1:N;
t2=discreto2'/fs;
g2=A*sin(2*pi*frecuencia2*t2+fase);
g=[g1;g2];
t=[t1;t2];

d=ECG_signal1'+g+SigmaRuido*randn;
ref1=B*sin(2*pi*frecuencia1*t1);
ref2=B*sin(2*pi*frecuencia2*t2);
ref=[ref1;ref2];   
w = zeros(Nfiltro,N);

for k = Nfiltro:length(ECG_signal1)
    uk = flip(ref(k-Nfiltro+1:k,1));
    SenalEstimada(k)=uk'*w(:,k-1);
    w(:,k) = w(:,k-1) + mu * uk * (d(k) - SenalEstimada(k));
end

SenalLimpia=d-SenalEstimada;
error=ECG_signal1'-SenalLimpia;

figure (1)
plot(t,g,'m--','LineWidth',1.25);
hold on;
plot(t,SenalEstimada);
hold on;
plot(t,ref);
grid on;
legend('Tono puro','Tono estimado','Referencia');
title('Tonos f_T=60/70 Hz');
axis([46.525 46.575 -1.15 1.55]);
% print('ejercicio4b_tonos.png','-dpng');


figure (2)
plot(t,ECG_signal1);
hold on;
plot(t,SenalLimpia);
grid on;
legend('Señal ECG','Señal ECG estimada');
title('Señales de ECG real y estimada');
xlim([45.5 49]);
% print('ejercicio4b_ecg.png','-dpng');

figure (3)
plot(t,error);
grid on;
title('Error entre las señales real y estimada');
axis([0 80 -0.6 0.6]);
% print('ejercicio4b_error.png','-dpng');