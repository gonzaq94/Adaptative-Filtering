clear all;
close all;
load senal_ECG1.mat
format long

A=1;
B=0.1*A;
SNR=20;
SigmaRuido= sqrt (A^2 / (2*10^(SNR/10)));
fs=8192;
frecuencia=60;
N=length(ECG_signal1);
fase=pi/2;
Nfiltro=100;
mu=0.001;

discreto=1:N;
t=discreto'/fs;
g=A*sin(2*pi*frecuencia*t+fase);

d=ECG_signal1'+g+SigmaRuido*randn;
ref=B*sin(2*pi*frecuencia*t);
   
w = zeros(Nfiltro,N);

for k = Nfiltro:length(ECG_signal1)
    uk = flip(ref(k-Nfiltro+1:k,1));
    w(:,k) = w(:,k-1) + mu * uk * (d(k) - uk'*w(:,k-1));
end

SenalEstimada=filter(w(:,N),1,ref);
SenalLimpia=d-SenalEstimada;
error=ECG_signal1'-SenalLimpia;

figure (1)
plot(t,g,'r--','LineWidth',1.5);
hold on;
plot(t,SenalEstimada);  
hold on;
plot(t,ref);
grid on;
legend('Tono puro','Tono estimado','Referencia');
title('Tonos f_T=60 Hz');%Tonos puro, de referencia y estimado
axis([17.45 17.5 -1.15 1.5]);
% print('ejercicio4a_tonos.png','-dpng');

figure (2)
plot(t,ECG_signal1);
hold on;
plot(t,SenalLimpia);
grid on;
legend('Señal ECG','Señal ECG estimada');
title('Señales de ECG real y estimada');
% xlim([17 18]);
% print('ejercicio4a_ecg.png','-dpng');

figure (3)
plot(t,error);
grid on;
axis([0 0.12 -0.14 0.25]);
title('Error entre las señales real y estimada');
% print('ejercicio4a_error.png','-dpng');