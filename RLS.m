close all;
clear all;
load ir_short.mat
format long

SNR=5;
sigma_ruido = 1;
a = 0.9;
M = length(w0);
kvec = 0:(M-1);
r = (sigma_ruido^2/(1 - a^2)) * a.^abs(kvec);
Ru = toeplitz(r);
aval_Ru = eig(Ru);
sigma_v = sqrt(w0'*Ru*w0/10^(SNR/10));

n_it = 10;
u0 = 0;
N = 800;
rho=0.0001;
LambdaVectorRLS=[1, 0.99,0.98];
mismatchRLS=zeros(3,N);
M=length(w0);


W0 = zeros(size(w0));
for i = 1:N
    W0(:,i) = w0;
end

for itK=1:length(LambdaVectorRLS)
    
    Pk=1/rho * eye(M,M);
    lambda=LambdaVectorRLS(itK);

    deltaRLS = zeros(1,N);

    for i = 1:n_it

        u = zeros(N,1);
        u(1) = u0;
        for m = 2:N
            u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
        end

        M = length(w0);
        w = zeros(M,N);
        wRLS = zeros(M,N);
        y = filter(w0,1,u);

        for k = M:length(u)
            uk = flip(u(k-M+1:k,1));
            d=y(k)+sigma_v*randn;
            Pk=1/lambda * ( Pk - 1/lambda * Pk * ( uk * uk' ) * Pk / ( 1 +  1/lambda * uk' * Pk * uk) );
            wRLS(:,k) = wRLS(:,k-1) + Pk * uk * (d - uk' * wRLS(:,k-1));
        end

        deltaRLS = deltaRLS + sum((W0-wRLS).^2);
        end

    mismatchRLS(itK,:) = 10*log10(deltaRLS/(n_it*w0'*w0));

end

figure (1)
for i=1:length(LambdaVectorRLS)
    plot(mismatchRLS(i,:));
    hold on;
end
title('Algoritmos LMS,APA y RLS con a=0.9 y \mu optimizado');
grid on;
xlabel('tiempo discreto(n)');
ylabel('mismatch');
legend('RLS \lambda=0.89','RLS \lambda=0.9','RLS \lambda=0.91','Location','Southeast');
