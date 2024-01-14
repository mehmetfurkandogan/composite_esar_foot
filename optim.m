%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;

core_opt = false;

%%
laminaMin = 10;
laminaMax = 30;

Tbest = cell(1, laminaMax - laminaMin + 1);
SRbest = zeros(1, laminaMax - laminaMin + 1);
exitflag = zeros(1, laminaMax - laminaMin + 1);

opts = optimoptions(@ga, ...
                    'PopulationSize', 100, ...
                    'MaxGenerations', 200, ...
                    'EliteCount', 5, ...
                    'FunctionTolerance', 1e-2, ...
                    'PlotFcn', @gaplotbestf,...
                    'Display','iter');

theta = cell(1,laminaMax-laminaMin + 1);
for iter = 1:laminaMax-laminaMin + 1
    %% Constraints
    laminaCount = iter + laminaMin - 1; % 2 additional parameters for core thickness and material
    if core_opt == true
        nvars = laminaCount + 2;
        lb = [-ones(1, laminaCount) 0 0];
        ub = [2*ones(1, laminaCount) 1 3];
    else
        nvars = laminaCount + 1;
        lb = [-ones(1, laminaCount) 0];
        ub = [2*ones(1, laminaCount) 1];
    end

    %% Optimization
    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        nvars, [], [], [], [], lb, ub, [], 1:nvars, opts);
    toc
    
    if core_opt == true
        core = Tbest{iter}(end);
        material_id = Tbest{iter}(end-1);
        theta_down = Tbest{iter}(1:stack)*45;
        theta_up = Tbest{iter}(stack+1:end-2)*45;
        if core ~= 0
            load("Materials\Rohacell.mat")
        end
    else
        core = 0;
        material_id = Tbest{iter}(end);
        theta_down = Tbest{iter}(1:stack)*45;
        theta_up = Tbest{iter}(stack+1:end-1)*45;
    end
    if material_id == 1 
        material = "Cycom 381 IM7 UD";
        load("Materials\Cycom 381 IM7 UD.mat")
    else 
        material = "Carbone TWILL 200 gsm";
        load("Materials\Carbone TWILL 200 gsm.mat")
    end
    stack = round(laminaCount/3);

    %% Mass of the foot
    % Foot Dimensions
    shoe_size = 42;     % eu
    L_data = 230e-3;    % m % Total Lenght of the foot
    L_model = ((shoe_size - 2 ) * 20 / 3)*1e-3;
    a = 0.31 * L_model; % Width of the foot
    n = laminaCount*2;

    area = L_model * a;
    if core ~= 0
        mass = rho * area * n*t + rho_core * area * 3 * core;
    else
        mass = rho * area * n*t;
    end
    
    %% Result
    fprintf('UPPER KNEEL: Optimum layer orientation for %d laminates:\t [%s]s \n',(laminaCount-stack)*2,join(string(theta_up), ','));
    fprintf('LOWER KNEEL: Optimum layer orientation for %d laminates:\t [%s]s \n',stack*2,join(string(theta_down), ','));
    fprintf('Material: %s\n', material);
    if core_opt == true
        fprintf('Core Thickness: %d mm \n', core);
        fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount*2,mass/SRbest(iter));
        theta{iter} = {laminaCount,mass/SRbest(iter),mass,Tbest{iter}(1:end-2)*45,material,Tbest{iter}(end)};
    else
        fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount*2,1/SRbest(iter));
        theta{iter} = {laminaCount,1/SRbest(iter),mass,Tbest{iter}(1:end-1)*45,material};
    end
end
save('theta.mat',"theta");