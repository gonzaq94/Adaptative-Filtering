close all;
clear all;
load ir_short.mat
format long

%Valor asintótico mismatch promedio para APA K=1,2,4,8

n_it = 5;
u0 = 0;
sigma_ruido = 1;
a = 0.95;
rho=0.0001;
M = length(w0);
mu_vector=[0.1 0.3 0.5 1 1.25 1.5  1.75 1.99 2];
Kvector=[1 2 4 8];
mismatch=zeros(length(Kvector),length(mu_vector));
MatrizIteraciones=zeros(length(Kvector),length(mu_vector));
variacion=1000*ones(length(Kvector),length(mu_vector));
SNR=20;
kvec = 0:(M-1);
r = (sigma_ruido^2/(1 - a^2)) * a.^abs(kvec);
Ru = toeplitz(r);
aval_Ru = eig(Ru);
sigma_v = sqrt(w0'*Ru*w0/10^(SNR/10));
VariacionMinima=n_it*sum((w0).^2)*sigma_v;

for itK=1:length(Kvector)
    K=Kvector(itK);
    if K~=1
        Uvector=zeros(K,M);
        Dvector=zeros(K,1);
    end
    
    for q=1:length(mu_vector)

        mu=mu_vector(q);
        for i = 1:n_it
            if mu==2
                Nmaximo=400;
            else
                Nmaximo=5*10^5;
            end
            delta=0;
            u = zeros(M,1);
            u(1) = u0;
            x = normrnd(0, sigma_ruido, [Nmaximo 1]);
            u = filter([1], [1, -a], x);
            w = zeros(M,1);
            y = filter(w0,1,u);
            
            k=M;
     
            while ( k<300 || variacion(itK,q)>VariacionMinima ) && k<Nmaximo
              
                
                uk = flip(u(k-M+1:k,1));
                d(k)=y(k)+sigma_v*randn;
                if K==1
                    Uvector=uk';
                    Dvector=d(k);
                else
                    Uvector=[uk'; Uvector(1:end-1,:)];
                    Dvector=[d(k); Dvector(1:end-1,:)];
                end
                w = w + mu * Uvector' /(rho*eye(K,K)+Uvector*Uvector')* (Dvector - Uvector*w);
                variacion(itK,q)=sum((w0-w).^2);
                k=k+1;
            end
            MatrizIteraciones(itK,q)=k;
            delta = delta + variacion(itK,q);
        end
        mismatch(itK,q) = 10*log10(delta/(n_it*w0'*w0));

    end
end

figure (6)
for i=1:length(Kvector)
    plot(mu_vector,mismatch(i,:));
    hold on;
end
legend('K=1','K=2','K=4','K=8','Location','Northwest');
title('Valor asintótico del mismatch para APA con a=0.95 SNR=20');
xlabel('\mu');
ylabel('Mismatch asintótico');
% print('ejercicio2APA_asintotico.png','-dpng');