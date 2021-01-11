function [Psnr] = mPSNR(D, D0)
sizD = size(D);
Psnr = zeros(1,sizD(1));
for i = 1:sizD(1)
    dim3 = sizD(2)*sizD(3);
    sub = squeeze(D(i,:,:,:));
    sub0 = squeeze(D0(i,:,:,:));
    sub = mat2tens(tens2mat(sub,3)',[32 32 dim3],3);
    sub0 = mat2tens(tens2mat(sub0,3)',[32 32 dim3],3);
    psnr = 0;
    for j = 1:dim3
        psnr = psnr + 10*log10(255^2/mse(sub(:, :, j) - sub0(:, :, j)));
    end
    Psnr(i) = psnr/dim3;
end
Psnr = mean(Psnr);
%reD = reshape(D,[sqrt(size(D,1)),sqrt(size(D,1)),size(D,2)*size(D,3)]);
%reD0 = reshape(D0,[sqrt(size(D,1)),sqrt(size(D,1)),size(D,2)*size(D,3)]);
%psnr = 0;%zeros(1,size(reD,3));
%for i = 1:size(reD,3)
%    psnr = psnr + 10*log10(255^2/mse(reD(:, :, i) - reD0(:, :, i)));
%end
%psnr = psnr/size(reD,3);
%subplot(1,2,1);
%imshow(reD0(:,:,500));
%subplot(1,2,2);
%imshow(reD(:,:,500));
end