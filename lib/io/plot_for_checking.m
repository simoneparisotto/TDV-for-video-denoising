stepstream = 10;

ii_slice = 5;
jj_slice = 5;
kk_slice = 3;

if true
    
    data_ij = squeeze(u(:,:,kk_slice));
    v1_ij   = squeeze(V(:,:,kk_slice,1));
    v2_ij   = squeeze(V(:,:,kk_slice,2));
    
    data_ik = squeeze(u(:,jj_slice,:));
    v3_ik   = squeeze(V(:,jj_slice,:,3));
    v4_ik   = squeeze(V(:,jj_slice,:,4));
    
    data_jk = squeeze(u(ii_slice,:,:));
    v5_jk   = squeeze(V(ii_slice,:,:,5));
    v6_jk   = squeeze(V(ii_slice,:,:,6));
    
    
    %% ATTEMPT
    if false

        %close(hs1)
        %close(hs2)
        
        stepstream_ij = 15;
        stepstream_k  = 30;
        
        hs0 = figure(1000);
        h = slice(uorig(:,:,:), [size(uorig,2)], [1], 1:stepstream_k:size(uorig,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_orig.png')
        
        hs1 = figure(1001);
        uplot = u(2:end-1,2:end-1,2:end-1);
        h = slice(uplot, [size(uplot,2)], [1], 1:stepstream_k:size(uplot,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_noise.png')
        
        hs2 = figure(1002);
        uplot = u(2:end-1,2:end-1,2:end-1);
        h = slice(uplot, [size(uplot,2)], [1], 1:stepstream_k:size(uplot,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        hold on
        % {i,j}
        iii = location.ii_v(1:stepstream_ij:end,1:stepstream_ij:end,1:stepstream_k:end);
        jjj = location.jj_v(1:stepstream_ij:end,1:stepstream_ij:end,1:stepstream_k:end);
        kkk = location.kk_v(1:stepstream_ij:end,1:stepstream_ij:end,1:stepstream_k:end);
        hhij = streamline(-V(:,:,:,1), V(:,:,:,2), zeros(size(uplot)-1),jjj,iii,kkk);
        set(hhij,'Color','blue');
        % {i,k}
        iii = location.ii_v(1:stepstream_ij:end,size(uplot,2)-2,1:stepstream_k:end);
        jjj = location.jj_v(1:stepstream_ij:end,size(uplot,2)-2,1:stepstream_k:end);
        kkk = location.kk_v(1:stepstream_ij:end,size(uplot,2)-2,1:stepstream_k:end);
        hhik = streamline(zeros(size(uplot)-1), V(:,:,:,4), -V(:,:,:,3), jjj,iii,kkk);
        set(hhik,'Color','red');
        % {j,k}
        iii = location.ii_v(2,1:stepstream_ij:end,1:stepstream_k:end);
        jjj = location.jj_v(2,1:stepstream_ij:end,1:stepstream_k:end);
        kkk = location.kk_v(2,1:stepstream_ij:end,1:stepstream_k:end);
        hhjk = streamline(V(:,:,:,6),zeros(size(uplot)-1), -V(:,:,:,5), jjj,iii,kkk);
        set(hhjk,'Color','yellow');
        %grid on; box on;
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_noise_with_v.png')
        
        hs4 = figure(1003);
        h = slice(aniso(:,:,:,1), [size(aniso,2)], [1], 1:stepstream_k:size(aniso,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_aniso_ij.png')
        
        hs5 = figure(1004);
        h = slice(aniso(:,:,:,2), [size(aniso,2)], [1], 1:stepstream_k:size(aniso,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_aniso_ik.png')
        
        hs6 = figure(1005);
        h = slice(aniso(:,:,:,3), [size(aniso,2)], [1], 1:stepstream_k:size(aniso,3));
        axis ij
        set(h, 'EdgeColor','none', 'FaceColor','interp')
        alpha(.6)
        colormap('gray')
        axis off
        view(-25, 8);
        set(gcf, 'Color', 'w');
        export_fig('./stack_aniso_jk.png')
        
    end
    
    %%
    % PLOT I-J
    [XXu,YYu] = meshgrid(1:1:N,M:-1:1);
    [XX,YY]   = meshgrid(1.5:1:(N-0.5),(M-0.5):-1:1.5);
  
    
    figure
    
    subplot(3,3,1)
    imagesc(XXu(:),YYu(:),data_ij)
    axis off
    axis image
    hold on
    streamline(XX,YY,-v2_ij,v1_ij,XX(1:stepstream:end,1:stepstream:end),YY(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('I-J gradient')
    
    subplot(3,3,4)
    imagesc(XXu(:),YYu(:),data_ij)
    axis off
    axis image
    hold on
    streamline(XX,YY,v1_ij,v2_ij,XX(1:stepstream:end,1:stepstream:end),YY(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('I-J gradient perp')
    
    subplot(3,3,7)
    imagesc(XX(:),YY(:),squeeze(aniso(:,:,kk_slice,1)))
    axis off
    axis image
    set(gca,'YDir','normal')
    title('anisotropy I-J')
    
    % PLOT I-K
    [TTu,XXu] = meshgrid(1:1:T,M:-1:1);
    [TT,XX]   = meshgrid(1.5:1:(T-0.5),(M-0.5):-1:1.5);
    
    subplot(3,3,2)
    imagesc(TTu(:),XXu(:),data_ik)
    axis off
    axis image
    hold on
    streamline(TT,XX,-v4_ik,v3_ik,TT(1:stepstream:end,1:stepstream:end),XX(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('I-K gradient')
    
    subplot(3,3,5)
    imagesc(TTu(:),XXu(:),data_ik)
    axis off
    axis image
    hold on
    streamline(TT,XX,v3_ik,v4_ik,TT(1:stepstream:end,1:stepstream:end),XX(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('I-K gradient perp')
    
    subplot(3,3,8)
    imagesc(TT(:),XX(:),squeeze(aniso(:,jj_slice,:,2)))
    axis off
    axis image
    set(gca,'YDir','normal')
    title('anisotropy I-K')
    
    % PLOT J-K
    [TTu,YYu] = meshgrid(1:1:T,N:-1:1);
    [TT,YY]   = meshgrid(1.5:1:T-0.5,(N-0.5):-1:1.5);
    
    subplot(3,3,3)
    imagesc(TTu(:),YYu(:),data_jk)
    axis off
    axis image
    hold on
    streamline(TT,YY,-v6_jk,v5_jk,TT(1:stepstream:end,1:stepstream:end),YY(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('J-K gradient')
    
    subplot(3,3,6)
    imagesc(TTu(:),YYu(:),data_jk)
    axis off
    axis image
    hold on
    streamline(TT,YY,v5_jk,v6_jk,TT(1:stepstream:end,1:stepstream:end),YY(1:stepstream:end,1:stepstream:end))
    set(gca,'YDir','normal')
    hold off
    title('J-K gradient perp')
    
    subplot(3,3,9)
    imagesc(TT(:),YY(:),squeeze(aniso(ii_slice,:,:,3)))
    axis off
    axis image
    set(gca,'YDir','normal')
    title('anisotropy J-K')
    
end