n_it = 1;
u0 = 0;
sigma_ruido = 1;
a = 0.95;
rho=0.0001;
VariacionMinima=n_it*sum((w0).^2)*0.000000000000000000000000000000000000000001;
Nmaximo=170000;
M = length(w0);
mu_vector=[0.2 0.5 1 1.25 1.5 1.7 1.8 1.9 2];
Kvector=[1 2 4 8];
mismatch=zeros(length(Kvector),length(mu_vector));
MatrizIteraciones=zeros(length(Kvector),length(mu_vector));
variacion=1000*ones(length(Kvector),length(mu_vector));
itK=1;
q=1;

    K=8;
    if K~=1
        Uvector=zeros(K,M);
        Dvector=zeros(K,1);
    end
    
 

        mu=1;
        for i = 1:n_it
            
            delta=0;
            u = zeros(M,1);
            u(1) = u0;
            for m = 2:Nmaximo
                u(m) = a*u(m-1) + sqrt(sigma_ruido)*randn;
            end

            w = zeros(M,1);
            y = filter(w0,1,u);
            
            k=M;
 
            while ( k<300 || variacion(itK,q)>VariacionMinima ) && k<Nmaximo
                
                uk = flip(u(k-M+1:k,1));
                if K==1
                    Uvector=uk';
                    Dvector=y(k);
                else
                    Uvector=[uk'; Uvector(1:end-1,:)];
                    Dvector=[y(k); Dvector(1:end-1,:)];
                end
                w = w + mu * Uvector' /(rho*eye(K,K)+Uvector*Uvector')* (Dvector - Uvector*w);
                variacion(itK,q)=sum((w0-w).^2);
                if k==800 
                    k=800;
                end
%                 u=[u ; a*u(k) + sqrt(sigma_ruido)*randn];
%                 y = filter(w0,1,u);
                k=k+1;
            end
            MatrizIteraciones(itK,q)=k;
            delta = delta + variacion(itK,q);
        end
        mismatch(itK,q) = 10*log10(delta/(n_it*w0'*w0));