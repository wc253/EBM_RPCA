% data: 30 subjects 11 view 21 illu imgsize 32*32
clc;
clear;close all;
addpath(genpath('lib'));
dataname = 'CMUface';
dataRoad = ['data/' dataname];
saveroad = ['result/result_for_' dataname];

%% Set enable bits
EN_HoRPCA       = 1;
EN_KBR_RPCA    = 1;
EN_TNN_RPCA   = 1;
EN_EBM_RPCA     = 1;  % the proposed method
methodname  = {'noisy','HoRPCA-S','KBR-RPCA','TNN-RPCA','EBM-RPCA'};
load(dataRoad);

%% experimental parameter settings with the variation of sparse noise ratio
sp_ratio = 0.2;                                                      % ratio of sparse noise
D0 = F(1:5,1:5,:,:);
sizeD = size(D0);
psnr = zeros(1,length(methodname));
deD = cell(1,length(methodname));
Time = zeros(1,length(methodname));

Par_tune = Parset(sp_ratio);

%% generate sparse salt&pepper noise
num = ceil(sp_ratio*prod(sizeD));
x = randperm(prod(sizeD),num);
D = D0;
D(x(1:num/2)) = 0;                                               % Minimum value
D(x((num/2+1):num)) = 1;                                    % Maximum  value
S = D - D0;

%% noise degree
enList = [];
j  = 1;
psnr(j) = mPSNR(255*D, 255*D0);
deD{j}  = D;
enList = [enList,j];

%% Use HoRPCA
j  = j + 1;
if EN_HoRPCA
    disp(['performing ',methodname{j}, ' ... ']);
    ParH.mu0         = 1/(length(sizeD)+1);
    ParH.mu1         = 5*std(D(:));
    ParH.mu2         = ParH.mu1;
    ParH.mu_min      = 1e-4;
    ParH.mu_max      = 1e2;
    ParH.max_iter    = 1000;
    ParH.lambdaS     = Par_tune{j};
    ParH.verbose     = false;
    ParH.V0          = cell(3,1);
    ParH.opt_tol     = 1e-5;
    for vi = 1:3
        ParH.V0{vi}  = zeros(32,32,sizeD(3));
    end
    ParH.X0          = tenzeros([32,32,sizeD(3)]);%B{1};
    ParH.E0          = tenzeros([32,32,sizeD(3)]);%F{1};
    ParH.lambda      = ParH.lambdaS./sqrt(max([32,32,sizeD(3)]));%params.lambdaS*0.25/sqrt(max(size(data.T)));
    psnr_c = [];
    deD_c = zeros(sizeD);
    tic;
    for s = 1:sizeD(1)
        for v = 1:sizeD(2)
            D0_c = reshape(squeeze(D0(s,v,:,:))',[32,32,sizeD(3)]);
            D_c = reshape(squeeze(D(s,v,:,:))',[32,32,sizeD(3)]);
            data.T           = D_c;
            data.X           = D0_c;
            ccc = tensor_rpca_adal(data, ParH);
            psnr_c(s,v) = cal_PSNR(255*double(ccc.X), 255*D0_c);
            deD_c(s,v,:,:) = double(tenmat(ccc.X,3));
        end
    end
    deD{j}       = deD_c;
    psnr(j) = mean(psnr_c(:));
    enList = [enList,j];
end

%% Use KBR-RPCA
j  = j + 1;
if EN_KBR_RPCA
    deD_c         = zeros(sizeD);
    psnr_c         = [];
    beta            = Par_tune{j}(1);%2.5*sqrt(max(sizeD));
    gamma       = Par_tune{j}(2);
    Par.maxIter = 1000;
    Par.lambda  = 1;%0.1;
    Par.mu        = 10;
    Par.tol         = 1e-5;
    Par.rho        = 1.05;
    disp(['performing ',methodname{j}, ' ... ']);
    tic;
    for s = 1:sizeD(1)
        for v = 1:sizeD(2)
            D_c = reshape(squeeze(D(s,v,:,:))',[32,32,sizeD(3)]);
            D0_c = reshape(squeeze(D0(s,v,:,:))',[32,32,sizeD(3)]);
            [ccc,~] = KBR_RPCA(D_c,beta,gamma,Par);
            psnr_c(s,v) = cal_PSNR(255*double(ccc), 255*D0_c);
            deD_c(s,v,:,:) = double(tenmat(ccc,3));
        end
    end
    deD{j} = deD_c;
    psnr(j) = mean(psnr_c(:));
    enList = [enList,j];
end

%% Use EN_TNN_RPCA
j  = j + 1;
if EN_TNN_RPCA
    deD_c = zeros(sizeD);
    psnr_c = [];
    TNN_lambda = Par_tune{j}*1/sqrt(sizeD(3)*32);
    opts.tol = 1e-8;
    opts.mu = 1e-4;
    opts.rho = 1.1;
    opts.DEBUG = 1;
    disp(['performing ',methodname{j}, ' ... ']);
    tic;
    for s = 1:sizeD(1)
        for v = 1:sizeD(2)
            D_c = reshape(squeeze(D(s,v,:,:))',[32,32,sizeD(3)]);
            D0_c = reshape(squeeze(D0(s,v,:,:))',[32,32,sizeD(3)]);
            [ccc, ~] = trpca_tnn(D_c,TNN_lambda,opts);
            psnr_c(s,v) = cal_PSNR(255*double(ccc), 255*D0_c);
            deD_c(s,v,:,:) = double(tenmat(ccc,3));
        end
    end
    deD{j} = deD_c;
    psnr(j) = mean(psnr_c(:));
    enList = [enList,j];
end

%% Use PBL-RPCA
j  = j + 1;
if EN_EBM_RPCA
    disp(['performing ',methodname{j}, ' ... ']);
    deD_c = zeros(sizeD);
    psnr_c = [];
    tic;
    for s = 1:sizeD(1)
        for v = 1:sizeD(2)
            D_c = reshape(squeeze(D(s,v,:,:))',[32,32,sizeD(3)]);
            D0_c = reshape(squeeze(D0(s,v,:,:))',[32,32,sizeD(3)]);
            S_c = reshape(squeeze(S(s,v,:,:))',[32,32,sizeD(3)]);
            [ccc, ~] = EBM_RPCA(D_c,100,D0_c,S_c,10^(-5),Par_tune{j}); 
            psnr_c(s,v) = cal_PSNR(255*double(ccc), 255*D0_c);
            deD_c(s,v,:,:) = double(tenmat(ccc,3));
        end
    end
    deD{j} = deD_c;
    psnr(j) = mean(psnr_c(:));
    enList = [enList,j];
end

%% Show result
%save('...');
show_FD(D0,deD,enList,methodname,psnr,sp_ratio)