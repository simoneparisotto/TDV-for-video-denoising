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
    % params = parameters for TDV
    [uorig,fps,M,N,C,T,I_MAX] = create_video(OBJECTS{ob});
    
    for nn = 1:numel(SIGMANOISES)
        
        % add noise to original video
        u = add_noise(uorig,SIGMANOISES(nn));
        % compute psnr u
        [PSNR_noise_final,PSNR_noise_fbf] = compute_psnr(u,uorig,I_MAX);
        
        for aa = 1:numel(ALPHAS)
            
            % CREATE DIR AND FILE
            result_dir = ['./results_with_linesearch_',date,'/',num2str(SIGMANOISES(nn)),'_',OBJECTS{ob},'_',num2str(ALPHAS{aa}(1)),'-',num2str(ALPHAS{aa}(2)),'-',num2str(ALPHAS{aa}(3)),'/'];
            if ~exist(result_dir,'dir')
                mkdir(result_dir);
            end
            fileID = fopen([result_dir,OBJECTS{ob},'_',num2str(SIGMANOISES(nn)),' ',datestr(now),'.txt'],'w');
            
            %% STARTING PARAMETERS
            ST_SIGMAS = 1;
            ST_RHOS   = 2;
            ETAS      = round(10*1./(20*(sigmanoise/255)^2))/10;
            
            PSNR_TDV_ETA = zeros(numel(sigmanoise),numel(ST_RHOS),numel(ETAS));
            verbose   = 1;
            maxiters  = 1000;
            tolerance = 1e-4;
            
            maxiter_linesearch = 500;
            experiments = [ST_SIGMAS, ST_RHOS, ETAS, NaN];
            step_s = 0.4;
            step_r = 0.4;
            step_e = 2;
            
            for ex = 1:size(experiments,1)
                params   = create_params_multiple(alpha,sigmanoise,fps,experiments(ex,1),experiments(ex,2),experiments(ex,3),maxiters,verbose,tolerance);
            end
            
            % params for run twice and improve
            params_rerun   = create_params(alpha,sigmanoise,fps,0.02,0.5,experiments(ex,3),maxiters,verbose,tolerance);
            
            fprintf(fileID,'\nObject:                  %s\n',object);
            fprintf(fileID,'Video size:              %dx%dx%dx%d\n',M,N,C,T);
            fprintf(fileID,'TDV lambda (order):      %s (%s)\n',num2str(params.lambda),num2str(params.order));
            fprintf(fileID,'Stdev:                   %2.4f%% (%d/255)\n',params.stdev,sigmanoise);
            fprintf(fileID,'Params iterations:       %03d\n',  params.niter);
            
            fprintf('\nObject:                  %s\n',object);
            fprintf('Video size:              %dx%dx%dx%d\n',M,N,C,T);
            fprintf('TDV lambda (order):      %s (%s)\n',num2str(params.lambda),num2str(params.order));
            fprintf('Stdev:                   %2.4f%% (%d/255)\n',params.stdev,sigmanoise);
            fprintf('Params iterations:       %03d\n\n',  params.niter);
            
            PSNR_TDV_fbf = zeros(size(experiments,1),T);
            
            u_TDV = zeros(M,N,C,T);
            
            for ex = 1:maxiter_linesearch
                
                params = create_params_multiple(alpha,sigmanoise,fps,experiments(end,1),experiments(end,2),experiments(end,3),maxiters,verbose,tolerance);
                
                [u_TDV,PSNR_TDV_fbf,PSNR_TDV_final,timecpu_TDV] = TDV_wrapper(u,uorig,params);
                experiments(end,4) = PSNR_TDV_final;
                
                stringstar = {'', '*'};
                fprintf(fileID,'TEST: %2d - (sigma,rho,eta) = (%2.2f,%2.2f,%2.2f) - PSNR TDV = %2.2f - cputime = %2.2f | %s\n',ex,params.sigma,params.rho,params.eta,PSNR_TDV_final,timecpu_TDV,stringstar{max(experiments(:,4)) == experiments(end,4)});
                fprintf('TEST: %2d - (sigma,rho,eta) = (%2.2f,%2.2f,%2.2f) - PSNR TDV = %2.2f - cputime = %2.2f | %s\n',ex,params.sigma,params.rho,params.eta,PSNR_TDV_final,timecpu_TDV,stringstar{max(experiments(:,4)) == experiments(end,4)});
                
                [experiments,flag_break] = new_fun_linesearch(experiments,step_s,step_r,step_e);
                
                tol_s = 0.05;
                tol_r = 0.05;
                tol_e = 0.25;
                
                while flag_break
                    
                    fprintf(fileID,'-- flag break recognized: all neighbours tested -- \n');
                    fprintf(fileID,'-- try to reduce the stepsize -- \n');
                    fprintf('-- flag break recognized: all neighbours tested -- \n');
                    fprintf('-- try to reduce the stepsize -- \n');
                    
                    % we reduce the timestep only if all the neighbourhoods
                    % are tested and there is no further direction to go
                    step_s = step_s/2;
                    step_r = step_r/2;
                    step_e = step_e/2;
                    
                    if (step_s<tol_s && step_r<tol_r && step_e<tol_e)
                        fprintf(fileID,'-- stepsizes below tolerance: stop the loop -- \n');
                        fprintf('-- stepsizes below tolerance: stop the loop -- \n');
                        flag_break = 1;
                        break
                    else
                        [experiments,flag_break] = fun_linesearch3(experiments,step_s,step_r,step_e);
                    end
                    
                end
                
                if flag_break
                    break
                end
                
            end
            
            % GET BEST PSNR and BUILD THE MATRIX
            [best_PSNR,ex_pos] = max(experiments(:,4));
            best_sigma         = experiments(ex_pos,1);
            best_rho           = experiments(ex_pos,2);
            best_eta           = experiments(ex_pos,3);
            
            if size(experiments,1)>1
                params = create_params(alpha,sigmanoise,fps,experiments(ex_pos,1),experiments(ex_pos,2),experiments(ex_pos,3),maxiters,verbose,tolerance);
                [u_TDV,PSNR_TDV_fbf,PSNR_TDV_final,timecpu_TDV] = TDV_wrapper(u,uorig,params);
                fprintf('BEST: %2d - (sigma,rho,eta) = (%2.2f,%2.2f,%2.2f) - PSNR TDV = %2.2f - cputime = %2.2f\n',ex_pos,experiments(ex_pos,1),experiments(ex_pos,2),experiments(ex_pos,3),experiments(ex_pos,4),timecpu_TDV);
            end
            
            %% ROF 2D+t
            [u_TV3D,PSNR_TV3D_fbf,PSNR_TV3D_final,timecpu_TV3D] = ROF2Dt_wrapper(u,uorig,I_MAX,params);
            fprintf('ROF 2D+t: eta = %2.2f - PSNR ROF 2+t = %2.2f - cputime = %2.2f\n',ETAS,PSNR_TV3D_final,timecpu_TV3D);      
       
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
            result_name = ['results_',num2str(sigmanoise),'_',object,'_',num2str(params.lambda(1)),'-',num2str(params.lambda(2)),'-',num2str(params.lambda(3)),'.mat'];
            save([result_dir,result_name],'object','sigmanoise','experiments','result_dir','u','uorig','u_TDV','timecpu_TDV','timecpu_VBM3D','timecpu_VBM4D','u_VBM3D','u_VBM4D','PSNR_noise_fbf','PSNR_VBM3D_fbf','PSNR_VBM4D_fbf','PSNR_TDV_fbf','PSNR_noise_final','PSNR_TDV_final','PSNR_VBM3D_final','PSNR_VBM4D_final','ST_SIGMAS','ST_RHOS','ETAS','params');
            
            %% PLOT FIGURES
            if plot_figures
                create_figures
                path_to_optimality
            end
            
            %% CLEAN
            close all
            clc
            clear PSNR_noise PSNR_TDV PSNR_VBM3D
            
            
        end
    end
end
