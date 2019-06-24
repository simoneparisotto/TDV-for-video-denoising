function [experiments,flag_break] = fun_linesearch(experiments,step_s,step_r,step_e)

flag_break = 0;

KN     = 21;
RADIUS = 2.5;
alpha  = 0.25;

increment = [0,0,0];

% SET OF THE NEIGHBOURHOODS
step(1,:)    = [step_s 0.00 0.00];
step(3,:)    = [0.00 step_r 0.00];
step(4,:)    = [0.00 0.00 step_e];
step(5,:)    = [-step_s 0.00 0.00];
step(6,:)    = [0.00 -step_r 0.00];
step(7,:)    = [0.00 0.00 -step_e];
step(8,:)    = [step_s step_r 0.00];
step(9,:)    = [step_s -step_r 0.00];
step(10,:)   = [-step_s +step_r 0.00];
step(11,:)   = [-step_s -step_r 0.00];
step(13,:)   = [step_s 0.00 step_e];
step(14,:)   = [step_s 0.00 -step_e];
step(15,:)   = [-step_s 0.00 step_e];
step(16,:)   = [-step_s 0.00 -step_e];
step(17,:)   = [0.00 step_r step_e];
step(18,:)   = [0.00 step_r -step_e];
step(19,:)   = [0.00 -step_r step_e];
step(20,:)   = [0.00 -step_r -step_e];

liststep = zeros(20,1);

% GET THE CURRENT MAX OF PSNR
[current_max_PSNR,current_max_pos] = max(experiments(:,4));
best_params                        = experiments(current_max_pos,1:3);

% IF AT LEAST 10 PSNR ARE ALREADY COMPUTED
if size(experiments,1)>10 && mod(size(experiments,1)+1,10)
    
    [idx, dist] = rangesearch(experiments(:,1:3),best_params,RADIUS);
    
    % delete itself from the results
    idx         = cell2mat(idx);  idx(1)  = [];
    dist        = cell2mat(dist); dist(1) = [];
    
    if ~isempty(idx)
        
        % get the nearest points
        neighb                  = experiments(idx,:);
        
        % compute all the gradients
        gradF                   = (current_max_PSNR - neighb(:,4))./dist.';
        
        % get the maximum of the gradients
        [max_grad, maxgrad_pos] = max(gradF(:));
        
        % get the increment
        if max_grad>0
            increment               = alpha*(best_params - neighb(maxgrad_pos,1:3))*max_grad;
        end
        
    end
    
else
    
    % CHOOSE ONE NEIGHBOUR RANDOMLY
    choosestep = randi(size(step,1),1);
    increment  = step(choosestep,:);
    
end

% get direction
direction               = round(100*(best_params + increment))/100;

%% CHEKC IF THE NEW DIRECTION IS ALREADY TESTED
tested = ismember(direction,experiments(:,1:3),'rows');

% if already tested  or the parameters violate the positivity try to perturbate and look at the neighbourhood
tol_pos = 1e-5;
while tested || direction(1)<=tol_pos || direction(2)<=tol_pos || direction(3)<=tol_pos || direction(1)>direction(2)
    
    % CHOOSE ONE NEIGHBOUR RANDOMLY
    choosestep = randi(size(step,1),1);
    if liststep(choosestep)==1
        % if the step is already choosen, select another one
        choosestep = randi(size(step,1),1);
    else
        % else flag the step tried
        liststep(choosestep) = 1;
    end
    direction  = best_params + step(choosestep,:);
    
    % check if direction is already tested
    tested = ismember(direction,experiments(:,1:3),'rows');
    
    % if all the neighbourhood are tested thenstop
    if isequal(liststep,ones(size(step,1),1))
        flag_break = 1;
        break
    end
    
end

% APPEND TO THE LIST
if ~flag_break
    experiments = [experiments;
        direction(1), direction(2), direction(3), NaN];
end

return