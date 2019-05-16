close all;
clear all;
load ir_short.mat;
format long;

n_it = 100; %numero de iteraciones de la simulacion de Montercarlo
mu_vector=[0.01;0.0095;0.009;0.0085;0.008;0.0075;0.007;0.0065;0.006;0.0055;0.005;0.0045;0.004;0.0035;0.003;0.0025;0.002;0.001;0.0005];
M = length(w0);
u0 = 0;
sigma_ruido = 1;
a = 0.3;
N = 20000;
SNR=20;
Misadjustment=zeros(length(mu_vector),1);
e=zeros(1,N); %Se utiliza para simular el Misadjustment
ei=zeros(1,N);

kvec = 0:(M-1);
r = (sigma_ruido^2/(1 - a^2)) * a.^abs(kvec);
Ru = toeplitz(r);
aval_Ru = eig(Ru);
sigma_v = sqrt(w0'*Ru*w0/10^(SNR/10));

Jmin=sigma_v^2;

W0 = zeros(size(w0));
for i = 1:N
    W0(:,i) = w0;
end
%WO contiene el wo óptimo repetido en todas sus columnas. Sirve para
%calcular la norma de (wo-w(k)) para todos los k.
    
for l=1:length(mu_vector)
    
    mu = mu_vector(l);
    delta = zeros(1,N);
    
    for i = 1:n_it
        
        u = zeros(N,1);
        u(1) = u0;
        for q = 2:N
            u(q) = a*u(q-1) + sigma_ruido*randn;
        end
        %Genero un vector de tamaño N que contiene el AR1
        w = zeros(M,N);
        y = filter(w0,1,u);
 
        for k = M:length(u)
            uk = flip(u(k-M+1:k,1));
            d=y(k)+sigma_v*randn;
            w(:,k) = w(:,k-1) + mu * uk * (d - uk'*w(:,k-1));
            ei(k)=(d - uk'*w(:,k-1))^2;
        end
        e=e+ei;
        delta = delta + sum((W0-w).^2);
    end
    e=e/n_it;
    e_corto=e(end-600:end);
    e_prom=sum(e_corto)/600;
    Misadjustment(l)=(e_prom-Jmin)/Jmin;
    e=zeros(size(e));
end



MisadjustmentTeorico=mu_vector/2*trace(Ru);

figure (2)
plot(mu_vector,Misadjustment);
hold on;
plot(mu_vector,MisadjustmentTeorico);
xlabel('\mu');
ylabel('misadjustment');
title('Algoritmo LMS con a=0.3 y SNR=20 dB');
legend('Simulado','Teórico');
grid on;
% print('ejercicio1Misadj.png','-dpng');