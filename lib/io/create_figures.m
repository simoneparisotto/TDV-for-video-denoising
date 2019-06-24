%% CREATE FIGURES

[M,N,C,T] = size(uorig);

% PARAMETERS
width  = 7;     % Width in inches
height = 7;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 25;      % Fontsize
lw     = 2.5;      % LineWidth
msz    = 15;       % MarkerSize
lcfsz = ceil(9*M/100);

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

% PSNR comparison
pause(1)
figure,
plot(1:size(u,4),PSNR_noise_fbf,'.-k',...
    1:size(u,4),PSNR_VBM3D_fbf,'.r-',...
    1:size(u,4),PSNR_VBM4D_fbf,'.m-',...
    1:size(u,4),PSNR_TDV_fbf(1,:),'.b-',...
    1:size(u,4),PSNR_TV3D_fbf(1,:),'.g-')
legend({'noisy','VBM3D','VBM4D','TDV','ROF 2D+t'},'position',[0.2 0.23 0.1 0.1],'fontsize',fsz)
set(gcf, 'Color', 'w');
ylabel('PSNR','Fontsize',fsz)
xlabel('Frames','Fontsize',fsz)
set(gca,'fontsize',fsz)
xlim([1,size(u,4)])
ylim([min(PSNR_noise_fbf(:))-1 max(max([reshape(PSNR_TV3D_fbf(1,:),[],1),PSNR_VBM3D_fbf(:),PSNR_VBM4D_fbf(:),PSNR_TDV_fbf(:)]))+1])
%title([OBJECTS{ob},' video - PSNR: TV3D > VBM3D in ',num2str(100*sum(PSNR_TV3D>PSNR_VBM3D)/numel(PSNR_TV3D)), '% of frames'],'Fontsize',fsz)
export_fig([result_dir,'PSNR_comparison_',OBJECTS{ob},'.png'])

% frame comparison
pause(1)
stack_orig  = [];
stack_noise = [];
stack_VBM3D = [];
stack_VBM4D = [];
stack_TDV   = [];
stack_TV3D  = [];
stack_frames = [5,10,15,20,25,30];
stack_frames = [40, 60 80,100,120];

for st = stack_frames
    stack_orig  = cat(2,stack_orig,uorig(:,:,1:C,st));
    stack_noise = cat(2,stack_noise,u(:,:,1:C,st));
    stack_VBM3D = cat(2,stack_VBM3D,u_VBM3D(:,:,1:C,st));
    stack_VBM4D = cat(2,stack_VBM4D,u_VBM4D(:,:,1:C,st));
    stack_TDV   = cat(2,stack_TDV,u_TDV(:,:,1:C,st));
    stack_TV3D  = cat(2,stack_TV3D,u_TV3D(:,:,1:C,st));
end

hstack_orig = figure;
imshow(stack_orig,[0,1]);
text(3,M-lcfsz-3,'ORIGINAL','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_orig.png'])

hstack_noise = figure;
imshow(stack_noise,[0,1]);
text(3,M-lcfsz-3,'NOISE','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_noise.png'])

hstackVBM3D = figure;
imshow(stack_VBM3D,[0,1]);
text(3,M-lcfsz-3,'VBM3D','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_VBM3D.png'])

hstackVBM4D = figure;
imshow(stack_VBM4D,[0,1]);
text(3,M-lcfsz-3,'VBM4D','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_VBM4D.png'])

hstack_TDV = figure;
imshow(stack_TDV,[0,1]);
text(3,M-lcfsz-3,'TDV','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_TDV.png'])

hstack_TV3D = figure;
imshow(stack_TV3D,[0,1]);
text(3,M-lcfsz-3,'ROF 2D+t','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_TV3D.png'])

hstack_all = figure;
lcfsz = ceil(9*M/100);
imshow(cat(1,stack_orig,stack_noise,stack_VBM3D,stack_VBM4D,stack_TDV,stack_TV3D),[0,1]);
text(3,1*M-lcfsz-3,'ORIGINAL','Color','green','FontSize',lcfsz,'FontWeight','bold')
text(3,2*M-lcfsz-3,'NOISE','Color','green','FontSize',lcfsz,'FontWeight','bold')
text(3,3*M-lcfsz-3,'VBM3D','Color','green','FontSize',lcfsz,'FontWeight','bold')
text(3,4*M-lcfsz-3,'VBM4D','Color','green','FontSize',lcfsz,'FontWeight','bold')
text(3,5*M-lcfsz-3,'TDV','Color','green','FontSize',lcfsz,'FontWeight','bold')
text(3,6*M-lcfsz-3,'ROF 2D+t','Color','green','FontSize',lcfsz,'FontWeight','bold')
export_fig([result_dir,'u_all.png'])


% OUTPUT VIDEO
if ~exist('plot_figures','var'), plot_figures = 1; end
if plot_figures
    u(u<0) = 0;
    u(u>1) = 1;
    u_VBM3D(u_VBM3D<0) = 0;
    u_VBM3D(u_VBM3D>1) = 1;
    u_VBM4D(u_VBM4D<0) = 0;
    u_VBM4D(u_VBM4D>1) = 1;
    u_TDV(u_TDV<0) = 0;
    u_TDV(u_TDV>1) = 1;
    u_TV3D(u_TV3D<0) = 0;
    u_TV3D(u_TV3D>1) = 1;
    writevideoMP4(uorig,  [result_dir,'/u_',OBJECTS{ob}],       params.FrameRate)
    writevideoMP4(u,      [result_dir,'/u_noise_',OBJECTS{ob}], params.FrameRate)
    writevideoMP4(u_VBM3D,[result_dir,'/u_VBM3D_',OBJECTS{ob}], params.FrameRate)
    writevideoMP4(u_VBM4D,[result_dir,'/u_VBM4D_',OBJECTS{ob}], params.FrameRate)
    writevideoMP4(u_TDV,  [result_dir,'/u_TDV_',OBJECTS{ob}],   params.FrameRate)
    writevideoMP4(u_TV3D, [result_dir,'/u_TV3D_',OBJECTS{ob}], params.FrameRate)
end

%             % PATH TO OPTIMALITY
%             [best_PSNR,ex_pos] = max(experiments(:,4));
%             step_ex_plot = 1:size(experiments,1)-1;
%             local_best   = zeros(numel(step_ex_plot),1);
%             current_best = -Inf;
%             for ex = step_ex_plot
%                 local_best_psnr = experiments(ex,4);
%                 if local_best_psnr > current_best
%                     local_best(ex) = 1;
%                     current_best = local_best_psnr;
%                 end
%             end
%             
%             
%             step_ex_plot = find(local_best).';
%             
%             %%
%             figure('Position',[0,0,1200,800]),
%             scatter3(experiments(:,1),experiments(:,2),experiments(:,3),200,experiments(:,4), 'fill')
%             grid on
%             hold on
%             plot3(experiments(step_ex_plot,1),experiments(step_ex_plot,2),experiments(step_ex_plot,3),'r.--')
%             xlabel('$\sigma$','Interpreter','latex')
%             ylabel('$\rho$','Interpreter','latex')
%             zlabel('$\eta$','Interpreter','latex')
%             set(gca,'fontsize',fsz)
%             xlim([ fix(10*(min(experiments(:,1))-0.1))/10, ceil(10*(max(experiments(:,1))+0.1))/10 ])
%             ylim([ fix(10*(min(experiments(:,2))-0.1))/10, ceil(10*(max(experiments(:,2))+0.1))/10 ])
%             zlim([ fix((min(experiments(:,3))-1))  , ceil((max(experiments(:,3))+1))   ])
%             set(gca,'xtick',fix(10*(min(experiments(:,1))-0.1))/10:0.2:ceil(10*(max(experiments(:,1))+0.1))/10)
%             set(gca,'ytick',fix(10*(min(experiments(:,2))-0.1))/10:0.2:ceil(10*(max(experiments(:,2))+0.1))/10)
%             set(gca,'ztick',fix((min(experiments(:,3))-1)):2:ceil((max(experiments(:,3))+1)))
%             set(gca,'XTickLabelRotation',50)
%             set(gca,'YTickLabelRotation',-30)
%             set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%             set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
%             set(gcf, 'Color', 'w');
%             
%             %set(gca,'fontsize',fsz)
%             %view([-65,25])
%             view([24,19])
%             
%             cl = colorbar;
%             set(cl, 'YTickLabel', cellstr(num2str(reshape(get(cl, 'YTick'),[],1),'%0.1f')) )
%             cl.Label.String = 'PSNR';
%             
%             axis square
%             box on
%             export_fig([result_dir,OBJECTS{ob},num2str(SIGMANOISES(nn)),'_TV3D_optimalitypath.png'])
%             