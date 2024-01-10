%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;
%%
laminaMin = 15;
laminaMax = 20;

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
    laminaCount = iter + laminaMin - 1;
    lb = -ones(1, laminaCount);
    ub = 2*ones(1, laminaCount);
    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        laminaCount, [], [], [], [], lb, ub, [], 1:laminaCount, opts);
    toc
    fprintf('Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(Tbest{iter}*45);
    fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount,1/SRbest(iter));
    theta{iter} = {laminaCount,1/SRbest(iter),Tbest{iter}*45};
end
save('theta.mat',"theta");