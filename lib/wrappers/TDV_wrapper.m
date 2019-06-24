function [u_TDV,PSNR_TDV_fbf,PSNR_TDV_final,timecpu] = TDV_wrapper(u,uorig,I_MAX,paramslocal)

%% TDV DENOISING

[M,N,C,T] = size(u);

u_TDV    = zeros(M,N,C,T);
timecpu  = zeros(C,1);
PSNR_TDV = zeros(C,T);

for cc=1:C
    
    uorig_local       = squeeze(uorig(:,:,cc,:));
    u_local           = squeeze(u(:,:,cc,:));
    
    [u_TDV_local,timecpu(cc)] = TDVtime_orders(u_local,uorig_local,paramslocal);
    u_TDV(:,:,cc,:)           = reshape(u_TDV_local,M,N,1,T);
    
    PSNR_TDV(cc,:)    = psnr3D(u_TDV_local,uorig_local);
    
end
timecpu        = sum(timecpu);

[PSNR_TDV_final,PSNR_TDV_fbf] = compute_psnr(u_TDV(:),uorig(:),I_MAX);

return

% addendum - try to avoid the boundary issue in time
%u     = u(:,:,:,[1 1 1:end end end]);
%uorig = uorig(:,:,:,[1 1 1:end end end]);

% addendum
%PSNR_TDV = PSNR_TDV(:,3:end-2);
    
% addendum - try to avoid the boundary issue in time
%uorig = uorig(:,:,:,[3:end-2]);
%u_TDV = u_TDV(:,:,:,[3:end-2]);
% I_MAX          = 1;
% PSNR_TDV_fbf   = sum(PSNR_TDV,1)/C;
% PSNR_TDV_final = 10*log10(I_MAX^2/mean((uorig(:)-u_TDV(:)).^2));
