close all;
clear all;
load ir_short.mat
format long
sigma_ruido = 1;
a = 0.9;
SNR=20;

M = length(w0);
kvec = 0:(M-1);
r = (sigma_ruido^2/(1 - a^2)) * a.^abs(kvec);
Ru = toeplitz(r);
aval_Ru = eig(Ru);
sigma_v = sqrt(w0'*Ru*w0/10^(SNR/10));

% %Algoritmos LMS
% n_it = 10;
% u0 = 0;
% N = 100000;
% mu = 0.0001;
% 
% W0 = zeros(size(w0));
% for i = 1:N
%     W0(:,i) = w0;
% end
% delta = zeros(1,N);
% for i = 1:n_it
%     
%     u = zeros(N,1);
%     u(1) = u0;
%     for m = 2:N
%         u(m) = a*u(m-1) + sigma_ruido*randn;
%     end
%     
%     M = length(w0);
%     w = zeros(M,N);
%     y = filter(w0,1,u);
%     
%     for k = M:length(u)
%         uk = flip(u(k-M+1:k,1));
%         d=y(k)+sigma_v*randn;
%         w(:,k) = w(:,k-1) + mu * uk * (d - uk'*w(:,k-1));
%     end
%     delta = delta + sum((W0-w).^2);
% end
% mismatch = 10*log10(delta/(n_it*w0'*w0));
%    
% figure (1)
% plot(mismatch');
% grid on;
% xlabel('tiempo discreto(n)');
% ylabel('mismatch');
% title('Algoritmo LMS con \mu=0.001, a=0.9 y SNR=20 dB');
% print('ejercicio2LMS.png','-dpng');

% %Algoritmo NLMS
% 
% n_it = 10;
% u0 = 0;
% N = 50000;
% rho=0.0001;
% mu_vector=[0.1 1];
% mismatch=zeros(length(mu_vector),N);
% 
% W0 = zeros(size(w0));
% for i = 1:N
%     W0(:,i) = w0;
% end
% 
% for q=1:length(mu_vector)
%     
%     mu=mu_vector(q);
%     delta = zeros(1,N);
%     for i = 1:n_it
% 
%         u = zeros(N,1);
%         u(1) = u0;
%         for m = 2:N
%             u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
%         end
% 
%         M = length(w0);
%         w = zeros(M,N);
%         y = filter(w0,1,u);
% 
%         for k = M:length(u)
%             uk = flip(u(k-M+1:k,1));
%             d=y(k)+sigma_v*randn;
%             w(:,k) = w(:,k-1) + mu/(rho+uk'*uk) * uk * (d - uk'*w(:,k-1));
%         end
%         delta = delta + sum((W0-w).^2);
%     end
%     mismatch(q,:) = 10*log10(delta/(n_it*w0'*w0));
%    
% end
% 
% figure (2)
% for i=1:length(mu_vector)
%     plot(mismatch(i,:));
%     hold on;
% end
% 
% grid on;
% xlabel('tiempo discreto (n)');
% ylabel('mismatch');
% title('Algoritmo NLMS con a=0.9 y SNR=20 dB');
% legend('\mu=0.1','\mu=1','Location','Southwest');
% print('ejercicio2NLMS.png','-dpng');

%Algoritmo APA

% n_it = 10;
% u0 = 0;
% sigma_ruido = 1;
% a = 0.95;
% N = 20000;
% rho=0.0001;
% mu_vector=[0.1 0.5 1];
% Kvector=[1 2 4];
% mismatch=zeros(length(Kvector),length(mu_vector),N);
% 
% W0 = zeros(size(w0));
% for i = 1:N
%     W0(:,i) = w0;
% end
% 
% for itK=1:length(Kvector)
%     K=Kvector(itK);
%     if K~=1
%         Uvector=zeros(K,M);
%         Dvector=zeros(K,1);
%     end
%     
%     for q=1:length(mu_vector)
% 
%         mu=mu_vector(q);
%         delta = zeros(1,N);
%         for i = 1:n_it
% 
%             u = zeros(N,1);
%             u(1) = u0;
%             for m = 2:N
%                 u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
%             end
% 
%             M = length(w0);
%             w = zeros(M,N);
%             y = filter(w0,1,u);
% 
%             for k = M:length(u)
%                 uk = flip(u(k-M+1:k,1));
%                 d=y(k)+sigma_v*randn;
%                 if K==1
%                     Uvector=uk';
%                     Dvector=d;
%                 else
%                     Uvector=[uk'; Uvector(1:end-1,:)];
%                     Dvector=[d; Dvector(1:end-1,:)];
%                 end
%                 w(:,k) = w(:,k-1) + mu * Uvector' *inv(rho*eye(K,K)+Uvector*Uvector')* (Dvector - Uvector*w(:,k-1));
%             end
%             delta = delta + sum((W0-w).^2);
%         end
%         mismatch(itK,q,:) = 10*log10(delta/(n_it*w0'*w0));
% 
%     end
% end
% 
% figure (3)
% m=zeros(length(mu_vector),N);
% for i=1:length(mu_vector)
%     m(i,:)=mismatch(1,i,:);
% end
% for i=1:length(mu_vector)
%     plot(m(i,:));
%     hold on;
% end
% grid on;
% xlabel('tiempo discreto (n)');
% ylabel('mismatch (dB)');
% title('Algoritmo APA con a=0.95, SNR=20 dB y K=1');
% legend('\mu=0.1','\mu=0.5','\mu=1','Location','Southwest');
% print('ejercicio2APA_K1.png','-dpng');
% 
% figure (4)
% m=zeros(length(mu_vector),N);
% for i=1:length(mu_vector)
%     m(i,:)=mismatch(2,i,:);
% end
% for i=1:length(mu_vector)
%     plot(m(i,:));
%     hold on;
% end
% grid on;
% xlabel('tiempo discreto (n)');
% ylabel('mismatch (dB)');
% title('Algoritmo APA con a=0.95, SNR=20 dB y K=2');
% legend('\mu=0.1','\mu=0.5','\mu=1','Location','Southwest');
% print('ejercicio2APA_K2.png','-dpng');
% 
% figure (5)
% m=zeros(length(mu_vector),N);
% for i=1:length(mu_vector)
%     m(i,:)=mismatch(3,i,:);
% end
% for i=1:length(mu_vector)
%     plot(m(i,:));
%     hold on;
% end
% grid on;
% xlabel('tiempo discreto (n)');
% ylabel('mismatch (dB)');
% title('Algoritmo APA con a=0.95, SNR=20 dB K=4');
% legend('\mu=0.1','\mu=0.5','\mu=1','Location','Southwest');
% print('ejercicio2APA_K4.png','-dpng');

%Valor asintótico mismatch promedio para APA K=1,2,4,8

% n_it = 1;
% u0 = 0;
% sigma_ruido = 1;
% a = 0.95;
% rho=0.0001;
% VariacionMinima=n_it*sum((w0).^2)*1e-30;
% Nmaximo=10e6;
% M = length(w0);
% mu_vector=[0.1 0.3 0.5 1 1.25 1.5  1.75 1.99 2];
% Kvector=[1 2 4 8];
% mismatch=zeros(length(Kvector),length(mu_vector));
% MatrizIteraciones=zeros(length(Kvector),length(mu_vector));
% variacion=1000*ones(length(Kvector),length(mu_vector));
% 
% for itK=1:length(Kvector)
%     K=Kvector(itK);
%     if K~=1
%         Uvector=zeros(K,M);
%         Dvector=zeros(K,1);
%     end
%     
%     for q=1:length(mu_vector)
% 
%         mu=mu_vector(q);
%         for i = 1:n_it
%             if mu==2
%                 Nmaximo=400;
%             else
%                 Nmaximo=10e6;
%             end
%             delta=0;
%             u = zeros(M,1);
%             u(1) = u0;
%             for m = 2:Nmaximo
%                 u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
%             end
% 
%             w = zeros(M,1);
%             y = filter(w0,1,u);
%             
%             k=M;
%      
%             while ( k<300 || variacion(itK,q)>VariacionMinima ) && k<Nmaximo
%                 
%                 uk = flip(u(k-M+1:k,1));
%                 if K==1
%                     Uvector=uk';
%                     Dvector=y(k);
%                 else
%                     Uvector=[uk'; Uvector(1:end-1,:)];
%                     Dvector=[y(k); Dvector(1:end-1,:)];
%                 end
%                 w = w + mu * Uvector' /(rho*eye(K,K)+Uvector*Uvector')* (Dvector - Uvector*w);
%                 variacion(itK,q)=sum((w0-w).^2);
%               
% %                 u=[u ; a*u(k) + sqrt(sigma_ruido)*randn];
% %                 y = filter(w0,1,u);
%                 k=k+1;
%             end
%             MatrizIteraciones(itK,q)=k;
%             delta = delta + variacion(itK,q);
%         end
%         mismatch(itK,q) = 10*log10(delta/(n_it*w0'*w0));
% 
%     end
% end
% 
% figure (6)
% for i=1:length(Kvector)
%     plot(mu_vector,mismatch(i,:));
%     hold on;
% end
% legend('K=1','K=2','K=4','K=8','Location','Northwest');
% title('Valor asintótico del mismatch para APA con a=0.95');
% xlabel('\mu');
% ylabel('Mismatch asintótico');
% print('ejercicio2APA_asintotico.png','-dpng');