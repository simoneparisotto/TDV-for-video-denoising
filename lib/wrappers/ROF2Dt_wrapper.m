function [u_TV3D,PSNR_TV3D_fbf,PSNR_TV3D_final,timecpu] = ROF2Dt_wrapper(u,uorig,I_MAX,paramslocal)

%% ROF 2D+t DENOISING

[M,N,C,T] = size(u);

u_TV3D    = zeros(M,N,C,T);
timecpu  = zeros(C,1);
PSNR_TV3D = zeros(C,T);

for cc=1:C
    
    uorig_local       = squeeze(uorig(:,:,cc,:));
    u_local           = squeeze(u(:,:,cc,:));
    
    [u_TV3D_local,timecpu(cc)] = ROF2Dt_orders(u_local,uorig_local,paramslocal);
    u_TV3D(:,:,cc,:)           = reshape(u_TV3D_local,M,N,1,T);
    
    PSNR_TV3D(cc,:)    = psnr3D(u_TV3D_local,uorig_local);
    
end
timecpu        = sum(timecpu);

[PSNR_TV3D_final,PSNR_TV3D_fbf] = compute_psnr(u_TV3D,uorig,I_MAX);
    
return

% addendum - try to avoid the boundary issue in time
%u     = u(:,:,:,[1 1 1:end end end]);
%uorig = uorig(:,:,:,[1 1 1:end end end]);

% addendum
%PSNR_TDV = PSNR_TDV(:,3:end-2);
% addendum - try to avoid the boundary issue in time
%uorig = uorig(:,:,:,[3:end-2]);
%u_TDV = u_TDV(:,:,:,[3:end-2]);
