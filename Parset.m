function Par = Parset(sp_ratio)
% we try to tune the parameters to make some comparison methods achieve the
% best performance
sp_ratio_arr = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6];   % different degrees of sparse noise
lamdaS_Ho = [1.6,1.2,1,0.8,0.8,0.6,0.6,0.6];
beta_KBR = [1e3,1e3,1e3,1e3,100,100,100,10];             
gamma_KBR = [1e5,1e5,1e5,1e5,1e4,1e4,1e4,1e2];
lambda_TNN = [2,2,1.6,1.2,1,1,0.8,0.8];          
alpha_PBL = [0.95,0.9,0.85,0.85,0.85,0.8,0.8,0.8];

[~,idx] = find(sp_ratio_arr==sp_ratio);
Par{2} = lamdaS_Ho(idx);
Par{3} = [beta_KBR(idx),gamma_KBR(idx)];
Par{4} = lambda_TNN(idx);
Par{5} = alpha_PBL(idx);
