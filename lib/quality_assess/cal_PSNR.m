%% Calculate mean psnr
function [Psnr] = cal_PSNR(D, D0)
sizD = size(D);
psnr = [];
count = 1;
for i =1:sizD(1)
    psnr(count) = 10*log10(255^2/mse(D(i,:, :) - D0(i,:, :)));
    count = count + 1;
end
Psnr = mean(psnr);
end