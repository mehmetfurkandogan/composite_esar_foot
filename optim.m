%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;
%%
laminaMin = 12;
laminaMax = 12;

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
    nvars = laminaCount + 2;
    lb = [-ones(1, laminaCount) 0 0];
    ub = [2*ones(1, laminaCount) 1 3];

    %% Optimization
    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        nvars, [], [], [], [], lb, ub, [], 1:nvars, opts);
    toc
    
    core = Tbest{iter}(end);
    if core ~= 0
        load("Materials\Rohacell.mat")
    end
    if Tbest{iter}(end-1) == 1 
        material = "Cycom 381 IM7 UD";
        load("Materials\Cycom 381 IM7 UD.mat")
    else 
        material = "Carbone TWILL 200 gsm";
        load("Materials\Carbone TWILL 200 gsm.mat")
    end
    stack = ceil(length(Tbest{iter}(1:end-2))/3);
    
    theta_up = Tbest{iter}(1:2*stack)*45;
    theta_down = Tbest{iter}(2*stack+1:end-2)*45;
    
    %% Mass of the foot
    % Foot Dimensions
    shoe_size = 42;     % eu
    L_data = 230e-3;    % m % Total Lenght of the foot
    L_model = ((shoe_size - 2 ) * 20 / 3)*1e-3;
    a = 0.31 * L_model; % Width of the foot
    n = laminaCount*2;

    area = L_model * a;
    if t_core ~= 0
        H = (n-2)*t+2*t_core;
        mass = rho * area * (H-2*t_core) + rho_core * area * 2 * t_core;
    else
        H = n*t;
        mass = rho * area * H;
    end
    
    %% Result
    fprintf('UPPER KNEEL: Optimum layer orientation for %d laminates:\t',2*stack);
    disp(theta_up);
    fprintf('LOWER KNEEL: Optimum layer orientation for %d laminates:\t',laminaCount-2*stack);
    disp(theta_down);
    fprintf('Material: %s\n', material);
    fprintf('Core Thickness: %d mm \n', Tbest{iter}(end));
    fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount,mass/SRbest(iter));
    theta{iter} = {laminaCount,mass/SRbest(iter),mass,Tbest{iter}(1:end-2)*45,Tbest{iter}(end-1),Tbest{iter}(end)};
end
save('theta.mat',"theta");