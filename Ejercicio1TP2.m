close all;
clear all;
load ir_short.mat;
format long;

n_it = 10; %numero de iteraciones de la simulacion de Montercarlo
mu_vector=[0.01;0.005;0.001;0.0005];
M = length(w0);
u0 = 0;
sigma_ruido = 1;
a = 0.3;
N = 20000;
SNR=20;
mismatch=zeros(length(mu_vector),N);
MismatchTeorico=zeros(length(mu_vector),N);
e=zeros(length(mu_vector),1); %Se utiliza para simular el Misadjustment

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
        end
        ei=(d - w(:,k)'*uk)^2; %Con el último w calculo la diferencia entre la salida real y la estimada para calcular
                                  % el Misadjustment
        e(l)=e(l)+ei;
        delta = delta + sum((W0-w).^2);
    end
    mismatch(l,:)= 10*log10(delta/(n_it*(norm(w0)^2)));
    
end

e=e/n_it; %Tomo el promedio
misadjustment=(e-Jmin)/Jmin;

for z=1:length(mu_vector)
    for i=1:N
        for l=1:M
            MismatchTeorico(z,i) = MismatchTeorico(z,i) + (mu_vector(z)*Jmin)/(2-mu_vector(z)*aval_Ru(l)) + (1-mu_vector(z)*aval_Ru(l))^(2*i)*...
            (abs(w0(l))^2 - (mu_vector(z)*Jmin)/(2-mu_vector(z)*aval_Ru(l)));
        end
    end
end

MismatchTeorico=10*log10(MismatchTeorico/(w0'*w0));

Colores=['b' 'r' 'g' 'm'];
figure (1)
for i=1:length(mu_vector)
    plot(mismatch(i,:),'Color',Colores(i),'LineStyle','--');
    hold on;
    plot(MismatchTeorico(i,:),'Color',Colores(i));
    hold on;
end

xlabel('tiempo discreto (n)');
ylabel('mismatch (dB)');
grid on;
legend('\mu=0.01 simulado','\mu=0.01 teórico','\mu=0.005 simulado','\mu=0.005 teórico','\mu=0.001 simulado',...
    '\mu=0.001 teórico','\mu=0.0005 simulado','\mu=0.0005 teórico','Location','Northeast');
title('Algoritmo LMS con a=0.3 y SNR=20 dB');
print('ejercicio1Mismatach.png','-dpng');

MisadjustmentTeorico=mu_vector/2*trace(Ru);

figure (2)
plot(mu_vector,misadjustment);
hold on;
plot(mu_vector,MisadjustmentTeorico);
xlabel('\mu');
ylabel('misadjustment');
title('Algoritmo LMS con a=0.3 y SNR=20 dB');
legend('Simulado','Teórico');
grid on;
print('ejercicio1Misadj.png','-dpng');