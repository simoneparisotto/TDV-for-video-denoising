function [PSNR_final,PSNR_fbf] = compute_psnr(u,uorig,I_MAX)

[M,N,C,T] = size(u);

PSNR = zeros(C,T);
for cc=1:C
    uorig_local  = squeeze(uorig(:,:,cc,:));
    u_local      = squeeze(u(:,:,cc,:));
    PSNR(cc,:)   = psnr3D(u_local,uorig_local);
end
PSNR_fbf   = sum(PSNR,1)/C;
PSNR_final = 10*log10(I_MAX^2/mean((uorig(:)-u(:)).^2));
