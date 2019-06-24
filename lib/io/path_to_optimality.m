% PATH TO OPTIMALITY
[best_PSNR,ex_pos] = max(experiments(:,4));
step_ex_plot = 1:size(experiments,1)-1;
local_best   = zeros(numel(step_ex_plot),1);
current_best = -Inf;
for ex = step_ex_plot
    local_best_psnr = experiments(ex,4);
    if local_best_psnr > current_best
        local_best(ex) = 1;
        current_best = local_best_psnr;
    end
end

step_ex_plot = find(local_best).';

figure('Position',[0,0,1200,800]),
scatter3(experiments(:,1),experiments(:,2),experiments(:,3),200,experiments(:,4), 'fill')
grid on
hold on
plot3(experiments(step_ex_plot,1),experiments(step_ex_plot,2),experiments(step_ex_plot,3),'r.--')
xlabel('$\sigma$','Interpreter','latex')
ylabel('$\rho$','Interpreter','latex')
zlabel('$\eta$','Interpreter','latex')
set(gca,'fontsize',fsz)
xlim([ fix(10*(min(experiments(:,1))-0.1))/10, ceil(10*(max(experiments(:,1))+0.1))/10 ])
ylim([ fix(10*(min(experiments(:,2))-0.1))/10, ceil(10*(max(experiments(:,2))+0.1))/10 ])
zlim([ fix((min(experiments(:,3))-1))  , ceil((max(experiments(:,3))+1))   ])
set(gca,'xtick',fix(10*(min(experiments(:,1))-0.1))/10:0.2:ceil(10*(max(experiments(:,1))+0.1))/10)
set(gca,'ytick',fix(10*(min(experiments(:,2))-0.1))/10:0.2:ceil(10*(max(experiments(:,2))+0.1))/10)
set(gca,'ztick',fix((min(experiments(:,3))-1)):2:ceil((max(experiments(:,3))+1)))
set(gca,'XTickLabelRotation',50)
set(gca,'YTickLabelRotation',-30)
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gcf, 'Color', 'w');

view([24,19])

cl = colorbar;
set(cl, 'YTickLabel', cellstr(num2str(reshape(get(cl, 'YTick'),[],1),'%0.1f')) )
cl.Label.String = 'PSNR';

for ex = step_ex_plot([1:9:end-9,end])
    lab =  sprintf('   %03d: (%02.2f,%02.2f,%02.2f)',ex,experiments(ex,1:3));
    txtplot = text(experiments(ex,1),experiments(ex,2),experiments(ex,3),lab,'Fontsize',18);
    uistack(txtplot, 'top')
end
axis square
box on
export_fig([result_dir,object,num2str(sigmanoise),'_TDV_optimalitypath.png'])
