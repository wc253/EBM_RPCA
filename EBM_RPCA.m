function [X_hat,lambda] = EBM_RPCA(Y,max_iters,X_ori,E_ori,std,alpha)
% Nonconvex Robust Low-Rank Tensor Reconstruction via an Empirical Bayes Method

% alpha needs to be tuned
N = size(Y);
n = prod(N);
% *** Control parameters ***
min_dX = 1*1e-6;
current_itr = 0;
lambda = std*10^5;%1e-1;%std+1e-4;
mode = length(N);

beta = 1;
for itr=1:mode
    Sigma{itr} = eye(N(itr));
    Sigma_inv{itr} = eye(N(itr));
end
Gamma = ones(prod(N),1);

max_SV = 2;
X_hat = zeros(n,1);
X_hat = tensor(reshape(X_hat,N));
E_hat = zeros(N);
z = zeros(n,2);

while (1)
    if(lambda>std)
        lambda = 0.8*lambda;
    end
  
    D_all = 1;
    for itr=1:mode
        [U{itr},D] = eig(Sigma{itr});
        UT{itr} = U{itr}';
        D = real(diag(D));
        D_all = kron(D,D_all);
    end

    z_old = z+1;
    loop_idx = 1;
    lambda0 = lambda*10^4;  % tune
    while((norm(z(:)-z_old(:),'fro')>1e-8)&&(loop_idx<200)) 
        
        if(lambda0>1*std)
            lambda0 = 0.9*lambda0;
        end
        L=2*max_SV/lambda0;
        D_all_inv = D_all./(2/L+ D_all);
        
        loop_idx = loop_idx+1;
        z_old = z;
        
        x_temp = z(:,1) -(2/lambda0/L)*sum(z,2) + (2/lambda0/L)*Y(:);
        e_temp = z(:,2) -(2/lambda0/L)*sum(z,2) + (2/lambda0/L)*Y(:);
        temp = full(ttensor(tensor(x_temp,N),UT));
        temp = temp.data;
        x = D_all_inv.*temp(:);
        X_hat = full(ttensor(tensor(reshape(x,N)),U));
        x = X_hat.data;
        x = x(:);
        
        e = Gamma./(Gamma+2/L).*e_temp;
        E_hat = reshape(e,N);
        
        z = [x e];
    end

% for itr00=1:10
    %% update phi
    for itr=1:mode
        [c,b] = eig((lambda/beta/2)^(1/mode)*eye(N(itr))+Sigma{itr});
%         [a,b,c] = svd((lambda)^(1/mode)/beta*eye(N(itr))+Sigma{itr});
        temp = find(diag(b)>1e-10);
        Grad{itr} = Sigma{itr} - Sigma{itr}*c(:,temp)*diag(1./(diag(b(temp,temp))))*c(:,temp)'*Sigma{itr};
        X_k{itr} = tenmat(X_hat,itr);
    end
    Grad_gamma = Gamma./(2*Gamma/lambda + 1);
    
    
    for itr = 1:mode%       
        Sigma_X_k = ttm(X_hat, Sigma_inv, setdiff(1:mode,itr));
        Sigma_X_k = tenmat(Sigma_X_k,itr);        
        
        Sigma{itr} = (alpha*N(itr)/prod(N))*double(Sigma_X_k)*double(X_k{itr})'+Grad{itr}; 
        [c,b] = eig(Sigma{itr});
        temp = find(diag(b)>1e-10);
        Sigma_inv{itr} = c(:,temp)*diag(1./(diag(b(temp,temp))))*c(:,temp)';
    end
    Gamma = (1-alpha)*E_hat(:).^2+Grad_gamma; 


    if(current_itr==0)
        dX = 1;        
    else
        temp = full(X_hat);
        dX = full(X_old)-temp;
        dX = norm(dX(:),'fro')/norm(temp(:),'fro');
%         error = X_ori-double(temp);
%         re = norm(error(:),'fro')/norm(X_ori(:),'fro');
%         error_E = E_ori-E_hat;
%         re_E = norm(error_E(:),'fro')/norm(E_ori(:),'fro');
%         disp([sprintf('Ier: %.1f error=%.2f error_E=%.2f.',current_itr+1,re,re_E)])
    end
    
    
    current_itr = current_itr+1;
    if (current_itr >= max_iters) 
        break  
    end
    if (dX < min_dX) 
        break
    end
    X_old = X_hat;
end