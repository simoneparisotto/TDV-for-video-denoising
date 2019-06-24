function [u_VBM3D,timecpu_VBM3D,PSNR_VBM3D_final,PSNR_VBM3D_fbf] = VBM3D_wrapper(u,uorig,sigma,I_MAX)

[M,N,C,T]      = size(u);
u_VBM3D        = zeros(M,N,C,T);
PSNR_VBM3D_fbf = zeros(C,T);

timecpu_VBM3D = cputime;
for cc=1:C
    uorig_local       = squeeze(uorig(:,:,cc,:));
    u_local           = squeeze(u(:,:,cc,:));
    
    [~, u_VBM3D_local] = VBM3D(u_local,sigma,T,0,uorig_local,'np');
    u_VBM3D(:,:,cc,:)  = double(reshape(u_VBM3D_local,M,N,1,T));
    
    % COMPARE PSNR frame-by-frame
    PSNR_VBM3D_fbf(cc,:)  = psnr3D(double(u_VBM3D_local),uorig_local);
    
    clear uorig_local u_local u_VBM3D_local
end
timecpu_VBM3D = cputime-timecpu_VBM3D;

[PSNR_VBM3D_final,PSNR_VBM3D_fbf] = compute_psnr(u_VBM3D,uorig,I_MAX);
        
return