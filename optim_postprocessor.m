%% Optimization Postprocesser
clc;clear;close all;
%%
load('theta.mat');
number_of_layers = zeros(size(theta));
strength_ratio = zeros(size(theta));
weight = zeros(size(theta));
for i = 1:length(theta)
    number_of_layers(i) = theta{i}{1};
    strength_ratio(i) = theta{i}{2};
    weight(i) = theta{i}{3};    % kg
end

figure;
hold on;grid on;
plot(number_of_layers,strength_ratio,LineWidth=1.5);
xticks(number_of_layers);
ylabel('Strength Ratio');
yyaxis right
plot(number_of_layers,weight*1e3,LineWidth=1.5);
ylabel('Weight (g)');
xlabel('Number of Layers');

