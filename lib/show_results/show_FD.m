function show_FD(D0,deD,enList,methodname,psnr,sp_ratio)
sizD = size(D0);
faces = 10;
ill = 10;
f_po = randperm(prod(sizD(1:2)),faces);
n = length(enList);
I0 = zeros(32,faces*32);
I_D = I0;
D0 = reshape(squeeze(D0(:,:,ill,:)),[sizD(1)*sizD(2),sizD(4)]);
D0 = reshape(D0(f_po,:)',[32 32 faces]);
plot_m  = tight_subplot(n+1,1,[.05 .03],[.01 .05],[.01 .01]);
for i = 1:faces
    I_D(:,(i*32-31):(i*32))           = D0(:,:,i);
end
axes(plot_m(1));
imshow(I_D);title('Original');
for f = 2:n+1
    I_D = I0;
    deDf = deD{enList(f-1)};
    D = reshape(squeeze(deDf(:,:,ill,:)),[sizD(1)*sizD(2),sizD(4)]); 
    D = reshape(D(f_po,:)',[32 32 faces]);
    for i = 1:faces
        I_D(:,(i*32-31):(i*32))           = D(:,:,i);
    end
    axes(plot_m(f));
    imshow(I_D);
    if f == 2
        title([methodname{enList(f-1)} ' (' num2str(100*sp_ratio) '% outlier) PSNR=' num2str(psnr(enList(f-1)),'%.2f') 'dB']);
    else
        title([methodname{enList(f-1)} ' PSNR=' num2str(psnr(enList(f-1)),'%.2f') 'dB']);
    end
end
end