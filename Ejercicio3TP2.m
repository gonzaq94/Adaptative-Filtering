close all;
clear all;
load ir_short.mat
format long

SNR=5;
sigma_ruido = 1;
a = 0;
M = length(w0);
kvec = 0:(M-1);
r = (sigma_ruido^2/(1 - a^2)) * a.^abs(kvec);
Ru = toeplitz(r);
aval_Ru = eig(Ru);
sigma_v = sqrt(w0'*Ru*w0/10^(SNR/10));

n_it = 10;
u0 = 0;
N = 20000;
rho=0.0001;
mu_vector=[1 1 1];
muLMS=0.0009;
lambda=0.98;
Kvector=[1 2 4];
ControlLMS=1;
mismatchAPA=zeros(length(Kvector),N);
mismatchLMS=zeros(1,N);
mismatchRLS=zeros(1,N);
M=length(w0);
Pk=1/rho * eye(M,M);

W0 = zeros(size(w0));
for i = 1:N
    W0(:,i) = w0;
end

for itK=1:length(Kvector)
    K=Kvector(itK);
    if K~=1
        Uvector=zeros(K,M);
        Dvector=zeros(K,1);
    end
    
    mu=mu_vector(itK);
    delta = zeros(1,N);
    deltaLMS = zeros(1,N);
    deltaRLS = zeros(1,N);

    for i = 1:n_it

        u = zeros(N,1);
        u(1) = u0;
        for m = 2:N
            u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
        end

        M = length(w0);
        w = zeros(M,N);
        wLMS = zeros(M,N);
        wRLS = zeros(M,N);
        y = filter(w0,1,u);

        for k = M:length(u)
            uk = flip(u(k-M+1:k,1));
            d=y(k)+sigma_v*randn;
            if K==1
                Uvector=uk';
                Dvector=d;
            else
                Uvector=[uk'; Uvector(1:end-1,:)];
                Dvector=[d; Dvector(1:end-1,:)];
            end
            w(:,k) = w(:,k-1) + mu * Uvector' *inv(rho*eye(K,K)+Uvector*Uvector')* (Dvector - Uvector*w(:,k-1));
            if ControlLMS==1
                wLMS(:,k) = wLMS(:,k-1) + muLMS * uk * (d - uk'*wLMS(:,k-1));            
                Pk=1/lambda * ( Pk - 1/lambda * Pk * ( uk * uk' ) * Pk / ( 1 +  1/lambda * uk' * Pk * uk) );
                wRLS(:,k) = wRLS(:,k-1) + Pk * uk * (d - uk' * wRLS(:,k-1));
            end
        end

        delta = delta + sum((W0-w).^2);
        
        if ControlLMS==1
            deltaLMS = deltaLMS + sum((W0-wLMS).^2);
            deltaRLS = deltaRLS + sum((W0-wRLS).^2);
        end
        end
    mismatchAPA(itK,:) = 10*log10(delta/(n_it*w0'*w0));
   
    if ControlLMS==1
        mismatchLMS(:) = 10*log10(deltaLMS/(n_it*w0'*w0));
        mismatchRLS(itK,:) = 10*log10(deltaRLS/(n_it*w0'*w0));
    end
    ControlLMS=0;
end

figure (1)
for i=1:length(Kvector)
    plot(mismatchAPA(i,:));
    hold on;
end
plot(mismatchRLS(:));
hold on;
plot(mismatchLMS(:));
title('Algoritmos LMS,APA y RLS con a=0, SNR=5 dB y \mu optimizado');
grid on;
xlabel('tiempo discreto (n)');
ylabel('mismatch (dB)');
legend('APA K=1 \mu=1','APA K=2 \mu=1','APA K=4 \mu=1','RLS \lambda=0.98','LMS \mu=0.0009','Location','Southeast');
print('ejercicio3b_blanco.png','-dpng');