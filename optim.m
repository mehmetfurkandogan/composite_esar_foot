%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;
%%
laminaMin = 9;
laminaMax = 9;

Tbest = cell(1, laminaMax - laminaMin + 1);
SRbest = zeros(1, laminaMax - laminaMin + 1);
exitflag = zeros(1, laminaMax - laminaMin + 1);

opts = optimoptions(@ga, ...
                    'PopulationSize', 10, ...
                    'MaxGenerations', 200, ...
                    'EliteCount', 5, ...
                    'FunctionTolerance', 1e-2, ...
                    'PlotFcn', @gaplotbestf,...
                    'Display','iter');

lb = [-ones(1, laminaCount) 0 0];
ub = [2*ones(1, laminaCount) 1 6];

theta = cell(1,laminaMax-laminaMin + 1);
for iter = 1:laminaMax-laminaMin + 1
    laminaCount = iter + laminaMin - 3; % 2 additional parameters for core thickness and material
    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        laminaCount, [], [], [], [], lb, ub, [], 1:laminaCount+2, opts);
    toc
    
    core = Tbest{iter}(end);
    material = Tbest{iter}(end-1);
    stack = ceil(length(Tbest{iter}(1:end-2))/3);
    
    theta_up = Tbest{iter}(1:2*stack)*45;
    theta_down = Tbest{iter}(2*stack+1:end-1)*45;

    fprintf('UPPER KNEEL: Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(theta_up);
    fprintf('LOWER KNEEL: Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(theta_down);
    fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount,1/SRbest(iter));
    theta{iter} = {laminaCount,1/SRbest(iter),Tbest{iter}(1:end-2)*45,Tbest{iter}(end-1),Tbest{iter}(end)};
end
save('theta.mat',"theta");