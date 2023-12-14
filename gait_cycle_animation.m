%% Gait Cycle Animation
% Mehmet Furkan Doğan
% 03 November 2023
clc;clear;close all;
%% IMPORT DATA
load('gait_cycle_data.mat','t','T','rib_x','rib_y','hip_x','hip_y','knee_x',...
    'knee_y','fibula_x','fibula_y','ankle_x','ankle_y','heel_x','heel_y',...
    'metat_x','metat_y','toe_x','toe_y',...
    'CoP','F_foot_ground_y','F_foot_ground_x','foot_theta');
% All positions are in mm
% Forces are in N
% Angles are in degrees
weight = 56.7*9.80665;  % Weight of the person in N
force_scale = 1/5;
ground_offset = 38; % mm
% ground_offset = 0; % mm
%% ANIMATION
% INITIALIZING THE FIGURE
f1 = figure('name','Gait Cycle','numberTitle','off');
hold on;
grid on;
set(gca,'NextPlot','replacechildren','DataAspectRatio',[1 1 1]);
xl = [1000 2000];
yl = [0 1200];
xlim(xl);ylim(yl);
xlabel('x (mm)');ylabel('y (mm)');
% f1.Position = [2         476        1361         208];
% INITIALIZING THE VIDEO
v1 = VideoWriter('Animations/walking_data.avi');
v1.FrameRate = 1/0.0145;
open(v1);
% DATA
figure(f1);
stance_phase_indices = 28:69;
spi = stance_phase_indices;
for i = stance_phase_indices
    HAT1_2_xr = [rib_x(i),hip_x(i)];
    HAT1_2_yr = [rib_y(i),hip_y(i)];
    Thigh_xr = [hip_x(i),knee_x(i)];
    Thigh_yr = [hip_y(i),knee_y(i)];
    Leg_xr = [fibula_x(i),ankle_x(i)];
    Leg_yr = [fibula_y(i),ankle_y(i)];
    Foot2_xr = [metat_x(i),toe_x(i)];
    Foot2_yr = [metat_y(i),toe_y(i)];
    Foot_xr = [heel_x(i),metat_x(i),ankle_x(i)];
    Foot_yr = [heel_y(i),metat_y(i),ankle_y(i)];
    foot = polyshape(Foot_xr,Foot_yr);
    plot(foot,LineWidth=1.5);
    hold on;
    plot(HAT1_2_xr,HAT1_2_yr,Thigh_xr,Thigh_yr,Leg_xr,Leg_yr,...
        Foot2_xr,Foot2_yr,'linewidth',1.5,'Color','k');
    ground = polyshape([xl(1)-100 xl(2)+100 xl(2)+100 xl(1)-100],...
                        [ground_offset,ground_offset,yl(1)-100,yl(1)-100]);
    plot(ground,LineWidth=1.5,FaceColor="#77AC30");
    quiver(CoP(i),ground_offset,F_foot_ground_x(i),F_foot_ground_y(i),force_scale,LineWidth=1.5,Color="#A2142F");
    movie_index = i-stance_phase_indices(1)+1;
    F1(movie_index) = getframe(gcf);
    writeVideo(v1,F1(movie_index));
    clf(f1)
    set(gca,'NextPlot','replacechildren','DataAspectRatio',[1 1 1]);
    xlabel('x (mm)');ylabel('y (mm)');
    xlim(xl);ylim(yl);
    grid on;
end
figure(f1);
% movie(F1,1,1/0.0145);
close(v1);
%% Plots
f2 = figure('name','Gait Cycle','numberTitle','off');
hold on;
grid on;
plot(100*t/T,100*F_foot_ground_x./weight,'r-',LineWidth=1.5);
plot(100*t/T,100*F_foot_ground_y./weight,'b-',LineWidth=1.5);
F_foot_ground = sqrt(F_foot_ground_x.^2 + F_foot_ground_y.^2);
plot(100*t/T,100*F_foot_ground./weight,'k-.',LineWidth=1.5);
legend('F_x','F_y','|F|')
xlabel('Percentage of the Gait Cycle (%)');
ylabel('Percentage of Total Weight (%)');
%% Coordinate Transformation
F_foot_ground_xp = zeros(106,1);
F_foot_ground_yp = zeros(106,1);
CoP_xp = zeros(106,1);
CoP_yp = zeros(106,1);

f3 = figure('name','Gait Cycle','numberTitle','off');
hold on;
grid on;
set(gca,'NextPlot','replacechildren','DataAspectRatio',[1 1 1]);
xl = [-50 250];
yl = [-50 150];
xlim(xl);ylim(yl);
xlabel('x'' (mm)');ylabel('y'' (mm)');
warning('off','MATLAB:polyshape:repairedBySimplify');
v3 = VideoWriter('Animations/foot_only.avi');
v3.FrameRate = 1/0.0145;
open(v3);
for i = stance_phase_indices
    theta = -atan2d(metat_y(i)-heel_y(i),metat_x(i)-heel_x(i));
    % theta = 90-atan2d(fibula_y(i)-ankle_y(i),fibula_x(i)-ankle_x(i));
    tx = heel_x(i);
    ty = heel_y(i);
    R = [cosd(theta) -sind(theta);
         sind(theta)  cosd(theta)];
    P = -R*[tx;ty];
    trans = [R(1,1) R(1,2) P(1);
             R(2,1) R(2,2) P(2);
             0      0      1 ];
    % heel = trans * [heel_x(i);heel_y(i);1];
    metat = trans * [metat_x(i);metat_y(i);1];
    heel = trans * [heel_x(i);heel_y(i);1];
    toe = trans * [toe_x(i);toe_y(i);1];
    ankle = trans * [ankle_x(i);ankle_y(i);1];
    fibula = trans * [fibula_x(i);fibula_y(i);1];
    CoP_t = trans * [CoP(i);ground_offset;1];
    CoP_xp(i) = CoP_t(1);
    CoP_yp(i) = CoP_t(2);
    ground1 = trans * [1000;ground_offset;1];
    ground2 = trans * [3000;ground_offset;1];
    Force = R * [F_foot_ground_x(i);F_foot_ground_y(i)];
    F_foot_ground_xp(i) = Force(1);
    F_foot_ground_yp(i) = Force(2);
    hold on;
    foot = polyshape([heel(1),metat(1),ankle(1)],[heel(2),metat(2),ankle(2)]);
    plot(foot,LineWidth=1.5);
    plot([metat(1),toe(1)],[metat(2),toe(2)],'k-',LineWidth=1.5);
    plot([ankle(1),fibula(1)],[ankle(2),fibula(2)],'k-',LineWidth=1.5);
    ground = polyshape(  [ground1(1) ground2(1) xl(2) xl(1)],...
                        [ground1(2) ground2(2) yl(1) yl(1)]);
    plot(ground,LineWidth=1.5,FaceColor="#77AC30");
    quiver(CoP_t(1),CoP_t(2),Force(1),Force(2),force_scale,LineWidth=1.5,Color="#A2142F");
    txt = compose('total: %.2f %%\nx: %.2f %%\ny: %.2f %%\n',...
                    100*sqrt(Force(1)^2+Force(2)^2)/weight,...
                    100*Force(1)/weight,100*Force(2)/weight);
    text(-40,125,txt)
    plot(CoP_xp(spi(1):i),CoP_yp(spi(1):i),'r-',LineWidth=1.5);
    movie_index = i-stance_phase_indices(1)+1;
    F3(movie_index) = getframe(gcf);
    writeVideo(v3,F3(movie_index));
    clf(f3)
    set(gca,'NextPlot','replacechildren','DataAspectRatio',[1 1 1]);
    xlim(xl);ylim(yl);
    xlabel('x'' (mm)');ylabel('y'' (mm)');
    xticks(xl(1):10:xl(end));
    grid on;
end
figure(f3);
% movie(F3,1,1/0.0145);
close(v3);
%% Plots
f4 = figure('name','Gait Cycle','numberTitle','off');
hold on;
grid on;
plot(100*t/T,...
    100*F_foot_ground_xp./weight,'r-',LineWidth=1.5);
plot(100*t/T,...
    100*F_foot_ground_yp./weight,'b-',LineWidth=1.5);
plot(100*t/T,...
    100*F_foot_ground./weight,'k-.',LineWidth=1.5);
legend('F_x''','F_y''','|F|')
xlabel('Percentage of the Gait Cycle (%)');
ylabel('Percentage of Total Weight (%)');
%% Center of Pressure
f5 = figure('name','Gait Cycle','numberTitle','off');
hold on;
grid on;
xlim([100*t(spi(1))/T,100*t(spi(end))/T]);
plot(100*t(spi)/T,CoP_xp(spi),'r-',LineWidth=1.5);
plot(100*t(spi)/T,CoP_yp(spi),'b-',LineWidth=1.5);
xlabel('Percentage of the Gait Cycle (%)');
ylabel('Center of Pressure Position (mm)');
legend('CoP_x''','CoP_y''',Location='best')
%%
f6 = figure('name','Gait Cycle','numberTitle','off');
set(gca,'NextPlot','replacechildren','DataAspectRatio',[1 1 1]);
hold on;
grid on;
plot(CoP_xp(spi),CoP_yp(spi),'r-',LineWidth=1.5);
scale = 0.01;
a = [CoP_xp(spi),CoP_yp(spi)]-[F_foot_ground_xp(spi),F_foot_ground_yp(spi)]*scale;
quiver(CoP_xp(spi),CoP_yp(spi),-F_foot_ground_xp(spi),-F_foot_ground_yp(spi),0.5,LineWidth=1.5,Color="#A2142F");
% plot(a(:,1),a(:,2),'r-',LineWidth=1.5);
xticks(xl(1):10:xl(end));
%% Save Data to a Table
gait_percentage = 100*t/T;
F_x_prime = 100*F_foot_ground_xp./weight;
F_y_prime = 100*F_foot_ground_yp./weight;
F_abs = 100*F_foot_ground./weight;
transformed_forces = table(gait_percentage,F_x_prime,F_y_prime,F_abs,...
    CoP_xp,CoP_yp);
writetable(transformed_forces,'transformed_forces.csv');