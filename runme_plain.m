%% TDV for Video Denoising
%  Copyright (c) 2019, Simone Parisotto
%  All rights reserved.
%
%  Author:
%  Simone Parisotto (email: sp751 at cam dot ac dot uk)
%      
%  Address:
%  Cambridge Image Analysis
%  Centre for Mathematical Sciences
%  Wilberforce Road
%  Cambridge CB3 0WA
%  United Kingdom
%  
%  Date:
%  March, 2019
%%

clear
close all
clc

delete(gcp('nocreate'))
warning off

addpath ./dataset
addpath ./lib/io
addpath ./lib/wrappers
addpath ./lib/lib_TDV/
addpath ./lib/BM3D_v2/
addpath ./lib/VBM4D_v1/
addpath ./lib/ROF2Dt/
addpath ./lib/export_fig-master

plot_figures  = 1;

%% CHOOSE VIDEO AND NOISE (in [0,255]
% select one or multiple videos among {'xylophone','tennis','flower','drop','coastguard','bus','foreman','salesman','function','xylophone_colour','drop_colour','function_colour'};
OBJECTS = {'drop'};
% select noise standard deviation (in the same intensity range of the
% video), e.g. one or multiple values among [10 20 35 50 70 90]
SIGMANOISES = [70];

% In case of multiple orders of derivatives (only [1 0 0] is currently implemented);
ALPHAS = {[1 0 0]};

for ob = 1:numel(OBJECTS)
    
    %% CREATE/LOAD ORIGINAL VIDEO
    % uorig  = real data
    % u      = noisy data
    [uorig,fps,M,N,C,T,I_MAX] = create_video(OBJECTS{ob});
    
    for nn = 1:numel(SIGMANOISES)
        
        % add noise to original video
        u = add_noise(uorig,SIGMANOISES(nn));
        
        % compute psnr u
        [PSNR_noise_final,PSNR_noise_fbf] = compute_psnr(u,uorig,I_MAX);
        
        for aa = 1:numel(ALPHAS)
            
            % CREATE DIR AND FILE
            result_dir = ['./results_without_linesearch_',date,'/',num2str(SIGMANOISES(nn)),'_',OBJECTS{ob},'_',num2str(ALPHAS{aa}(1)),'-',num2str(ALPHAS{aa}(2)),'-',num2str(ALPHAS{aa}(3)),'/'];
            if ~exist(result_dir,'dir')
                mkdir(result_dir);
            end
            fileID = fopen([result_dir,OBJECTS{ob},'_',num2str(SIGMANOISES(nn)),' ',datestr(now),'.txt'],'w');
            
            % STARTING PARAMETERS
            ETAS      = round(100*255/SIGMANOISES(nn))/100;
            ST_SIGMAS = 3.2/sqrt(ETAS);
            ST_RHOS   = 3.2/sqrt(ETAS);
            
            verbose   = 1;
            maxiters  = 1000;
            tolerance = 1e-3;
            
            experiments = [ST_SIGMAS, ST_RHOS, ETAS, NaN];
            
            params   = create_params(ALPHAS{aa},SIGMANOISES(nn),fps,experiments(1,1),experiments(1,2),experiments(1,3),maxiters,verbose,tolerance);
            
            fprintf(fileID,'\n');
            fprintf(fileID,'Object:                  %s\n',OBJECTS{ob});
            fprintf(fileID,'Video size:              %dx%dx%dx%d\n',M,N,C,T);
            fprintf(fileID,'TDV lambda (order):      %s (%s)\n',num2str(params.lambda),num2str(params.order));
            fprintf(fileID,'Stdev:                   %2.4f%% (%d/255)\n',params.stdev,SIGMANOISES(nn));
            fprintf(fileID,'Params iterations:       %03d\n',  params.niter);
            
            fprintf('\n');
            fprintf('Object:                  %s\n',OBJECTS{ob});
            fprintf('Video size:              %dx%dx%dx%d\n',M,N,C,T);
            fprintf('TDV lambda (order):      %s (%s)\n',num2str(params.lambda),num2str(params.order));
            fprintf('Stdev:                   %2.4f%% (%d/255)\n',params.stdev,SIGMANOISES(nn));
            fprintf('Params iterations:       %03d\n\n',  params.niter);
            
            %% TDV
            [u_TDV,PSNR_TDV_fbf,PSNR_TDV_final,timecpu_TDV] = TDV_wrapper(u,uorig,I_MAX,params);
            experiments(end,4) = PSNR_TDV_final;
            fprintf('TDV:      (sigma,rho,eta) = (%2.2f,%2.2f,%2.2f) - PSNR TDV = %2.2f - cputime = %2.2f\n',ST_SIGMAS,ST_RHOS,ETAS,PSNR_TDV_final,timecpu_TDV);
            
            %% ROF 2D+t
            [u_TV3D,PSNR_TV3D_fbf,PSNR_TV3D_final,timecpu_TV3D] = ROF2Dt_wrapper(u,uorig,I_MAX,params);
            fprintf('ROF 2D+t: eta = %2.2f - PSNR ROF 2+t = %2.2f - cputime = %2.2f\n',ETAS,PSNR_TV3D_final,timecpu_TV3D);      
            
        end
        
        %% V-BM3D
        % DISABLED: incompatible issue with R2018b+ APPLE OSX 10.14)
        u_VBM3D          = zeros(size(u));
        timecpu_VBM3D    = NaN;
        PSNR_VBM3D_final = NaN(size(PSNR_noise_final));
        PSNR_VBM3D_fbf   = NaN(size(PSNR_noise_fbf));
        % uncomment the line below for enabling V-BM3D:
        % [u_VBM3D,timecpu_VBM3D,PSNR_VBM3D_fbf,PSNR_VBM3D_final] = VBM3D_wrapper(u,uorig,SIGMANOISES(nn),I_MAX);
    
        %% V-BM4D
        [u_VBM4D,timecpu_VBM4D,PSNR_VBM4D_final,PSNR_VBM4D_fbf] = VBM4D_wrapper(u,uorig,SIGMANOISES(nn),I_MAX);
        
        %% PRINT RESULTS
        fprintf(fileID,'PSNR noise     : %2.2f\n',PSNR_noise_final);
        fprintf(fileID,'PSNR VBM3D     : %2.2f                     - cputime %2.2f s.\n',PSNR_VBM3D_final,timecpu_VBM3D);
        fprintf(fileID,'PSNR VBM4D     : %2.2f                     - cputime %2.2f s.\n',PSNR_VBM4D_final,timecpu_VBM4D);
        fprintf(fileID,'PSNR ROF 2D+t  : %2.2f (%2.2f,%2.2f,%2.2f) - cputime %2.2f s.\n',PSNR_TV3D_final,experiments(1,1),experiments(1,2),experiments(1,3),timecpu_TV3D);
        fclose(fileID);
        
        fprintf(       'PSNR noise:    %2.2f s.\n',PSNR_noise_final);
        fprintf(       'PSNR VBM3D:    %2.2f - cputime %2.2f s.\n',PSNR_VBM3D_final,timecpu_VBM3D);
        fprintf(       'PSNR VBM4D:    %2.2f - cputime %2.2f s.\n',PSNR_VBM4D_final,timecpu_VBM4D);
        fprintf(       'PSNR ROF 2D+t: %2.2f - cputime %2.2f s.\n',PSNR_TV3D_final,timecpu_TV3D);
        
        %% SAVE RESULTS
        pause(1)
        close all
        result_name = ['results_',num2str(SIGMANOISES(nn)),'_',OBJECTS{ob},'_',num2str(params.lambda(1)),'-',num2str(params.lambda(2)),'-',num2str(params.lambda(3)),'.mat'];
        save([result_dir,result_name]);
        
        %% PLOT FIGURES
        if plot_figures
            create_figures
        end
        
        %% CLEAN
        close all
        clc
        clear PSNR_noise PSNR_TDV PSNR_TV3D PSNR_VBM3D PSNR_VBM4D
        
    end
end