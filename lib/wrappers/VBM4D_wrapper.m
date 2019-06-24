function [u_VBM4D,timecpu_VBM4D,PSNR_VBM4D_final,PSNR_VBM4D_fbf] = VBM4D_wrapper(u,uorig,sigma,I_MAX)

[M,N,C,T]      = size(u);
u_VBM4D        = zeros(M,N,C,T);
PSNR_VBM4D_fbf = zeros(C,T);

timecpu_VBM4D = cputime;
for cc=1:C
    uorig_local       = squeeze(uorig(:,:,cc,:));
    u_local           = squeeze(u(:,:,cc,:));
    
    u_VBM4D_local      = vbm4d(u_local,sigma/255, 'np', 1, 1, 1, 0 );
    u_VBM4D(:,:,cc,:)  = double(reshape(u_VBM4D_local,M,N,1,T));
    
    % COMPARE PSNR frame-by-frame
    PSNR_VBM4D_fbf(cc,:)  = psnr3D(double(u_VBM4D_local),uorig_local);
    
    clear uorig_local u_local u_VBM4D_local
end
timecpu_VBM4D = cputime-timecpu_VBM4D;

[PSNR_VBM4D_final,PSNR_VBM4D_fbf] = compute_psnr(u_VBM4D,uorig,I_MAX);
    
return