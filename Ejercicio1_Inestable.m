close all;
clear all;
load ir_short.mat;
format long;

n_it = 20; %numero de iteraciones de la simulacion de Montercarlo
mu_vector=[0.015;0.014;0.013;0.012];
M = length(w0);
u0 = 0;
sigma_ruido = 1;
a = 0.3;
N = 20000;
SNR=20;
mismatch=zeros(length(mu_vector),N);

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
      
        delta = delta + sum((W0-w).^2);
    end
    mismatch(l,:)= 10*log10(delta/(n_it*(norm(w0)^2)));
   
end



figure (1)
for i=1:length(mu_vector)
    plot(mismatch(i,:));
    hold on;
end

xlabel('tiempo discreto (n)');
ylabel('mismatch (dB)');
grid on;
legend('\mu=0.015','\mu=0.014','\mu=0.013','\mu=0.012','Location','Northwest');
title('Divergencia del algoritmo LMS con a=0.3 y SNR=20 dB');
% print('ejercicio1inestable.png','-dpng');

display('La cota determinada por la small step size theory es:');
2/max(aval_Ru)

