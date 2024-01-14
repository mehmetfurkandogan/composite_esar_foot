%% CLT Calculations for an ESR Foot
% Mehmet Furkan Doğan
% 30.11.2023
clc;clear;close all;
%% Defining the material properties
core = false;

load('Materials/Cycom 381 IM7 UD.mat')

if core == true
    load("Materials\Rohacell.mat")
    t_core = 4e-3;
end
design_mass = 130;
shoe_size = 42;     % eu
upper = [0, 0, 0, 0, 0, 0, 0, -45, 45, 0, 0, 0, 90, 0, 0];
lower = [0, 0, 0, 0, 45, -45, 0, 0];
%% Calculation of compliance and stiffness matrices
% Compliance matrix for unidirectional lamina
S11 = 1/E_1;        % 1/Pa
S12 = -nu_12/E_1;   % 1/Pa
S13 = -nu_13/E_1;   % 1/Pa
S22 = 1/E_2;        % 1/Pa
S23 = -nu_23/E_2;   % 1/Pa
S33 = 1/E_3;        % 1/Pa
S44 = 1/G_23;       % 1/Pa
S55 = 1/G_31;       % 1/Pa
S66 = 1/G_12;       % 1/Pa
% Complete compliance matrix
S = [S11 S12 S13 0   0   0;
     S12 S22 S23 0   0   0;
     S13 S23 S33 0   0   0;
     0   0   0   S44 0   0;
     0   0   0   0   S55 0;
     0   0   0   0   0   S66];
% Stiffness matrix
C = inv(S);
% Reduced compliance matrix
S = [S11 S12 0;
     S12 S22 0;
     0   0   S66];  % 1/Pa
% Reduced stiffness matrix
Q = inv(S);         % Pa
%% Laminate Properties
tic
upper = [upper flip(upper)];
lower = [lower flip(lower)];
theta = [upper lower];

if t_core ~= 0
    theta_down = [theta(1:stack) 0 flip(theta(1:stack))];  % degree (Symmetric)
    theta_up = [theta(stack+1:end) 0 flip(theta(stack+1:end))];
    theta = [theta_up theta_down];
    n = size(theta,2);  % number of plies
    H = (n-2)*t+t_core_down+t_core_up;        % m % Total width of the lamimate
    h = zeros(1,n);
    h(1) = -H/2;
    core_up_id = length(theta_down) + 2*stack+1;
    core_down_id = stack + 1;
    for i = 1:length(theta)
        if i == core_up_id
            h(i+1) = h(i) + t_core_up;
        elseif i == core_down_id
            h(i+1) = h(i) + t_core_down;
        else
            h(i+1) = h(i) + t; % m
        end
    end
else
    % theta = [theta flip(theta)];                       % degree (Symmetric)
    n = size(theta,2);  % number of plies
    H = n*t;        % m % Total width of the lamimate
    for i = 0:length(theta)
        h(i+1) = -H/2 + i*t; % m
    end
    core_up_id = -1;
    core_down_id = -1;
end

%% Angle transformation
% Reuter matrix
R = [1 0 0;
     0 1 0;
     0 0 2];    % -

for i = 1:n
    if i == 1+ (n-1)/2 && core == true
        Qbar(:,:,i) = [1/E_core -nu_core/E_core 0;
                       -nu_core/E_core 1/E_core 0;
                       0 0 1/G_core];
        Sbar(:,:,i) = inv(Qbar(:,:,i));     % 1/Pa
        continue
    end
    s = sin(theta(i)*pi/180);
    c = cos(theta(i)*pi/180);
    % Transformation matrix
    T = [c^2    s^2     2*s*c;
         s^2    c^2     -2*s*c;
         -s*c   s*c     c^2-s^2];  % -
    Qbar(:,:,i) = inv(T)*Q*R*T*inv(R);  % Pa
    Sbar(:,:,i) = inv(Qbar(:,:,i));     % 1/Pa
end

%% Global laminate stiffness matrix
A = zeros(3,3); % Pa * m    % Extensional stiffnes matrix
B = zeros(3,3); % Pa * m^2  % Coupling stiffness matrix
D = zeros(3,3); % Pa * m^3  % Bending stiffness matrix

for i = 1:3
    for j = 1:3
        for k = 1:n
            A(i,j) = A(i,j) + Qbar(i,j,k) * (h(k+1)-h(k));
            B(i,j) = B(i,j) + (1/2) * Qbar(i,j,k) * (h(k+1)^2-h(k)^2);
            D(i,j) = D(i,j) + (1/3) * Qbar(i,j,k) * (h(k+1)^3-h(k)^3);
        end
    end
end

ABBD = [A B;
        B D];
%% Loading
% Foot Dimensions
L_data = 230e-3;    % m % Total Lenght of the foot
% L_model = 206e-3;  % m for proted design
L_model = ((shoe_size - 2 ) * 20 / 3)*1e-3;
b_rear = 0.311 * L_model;  % Lenght of the front part of the foot
b_front = L_model - b_rear;    % Lenght of the rear part of the foot
a = 0.31 * L_model; % Width of the foot
delta = 7.5;    % deg % the angle between foot axis and the walking direction
Iy = (1/12) * a*H^3; % m^4
Iz = (1/12) * H*a^3; % m^4
% J = Iy + Iz; % m^4
J = 1/16*a*H^3*(16/3-3.36*H/a*(1-H^4/(12*a^4)));

% Curved Beam
R = 100/208*L_model;
rn = H/log((R + H)/R);
e = R + H/2 - rn;
%% Total Mass of the Foot
area = L_model * a;
if core == true
    mass = rho * area * (H-t_core) + rho_core * area * t_core;
else
    mass = rho * area * H;
end

%% 
load('gait_forces.mat')

mass_data = 56.7;

number_of_time_steps = length(spi);
nots = number_of_time_steps;
% Toe contact
Fz = -F_foot_ground_yp(spi)*design_mass/mass_data;
Fy = F_foot_ground_xp(spi) * sind(delta)*design_mass/mass_data;
Fx = F_foot_ground_xp(spi) * cosd(delta)*design_mass/mass_data;
b = CoP_xp(spi)*1e-3 * L_model / L_data - b_rear;
% b = abs(b);
%% Loadings

Nx = Fx/(H*a)*H + Fy.*b*(a/2)*H/Iz + Fz.*b/(H*a*e)*(R*log((R+H/2+e)/(R-H/2+e)) - H);    % N/m
Ny = Fy./(H*b)*H;    % N/m
Nxy = Fy/(H*a)*H; %+ Fz*a^3/(8*J)*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2))));  % N/m
N = [Nx Ny Nxy]';

Mx = Fz.*b/(2*H*a*e)*(2*H*R + 2*R*(R+e)*log((R-H/2+e)/(R+H/2+e)));   % N*m/m
My = zeros(size(b));      % N*m/m
Mxy = Fz*a^4/(16*J)*(-1/4*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2)))) + 1/2*H/a*(H^2/a^2 + 1)^(3/2));                 % N*m/m
M = [Mx My Mxy]';

NM = [N;M];

%% Strains and curvatures
eps0kappa = ABBD\NM;

eps0 = eps0kappa(1:3,:);      % m/m
kappa = eps0kappa(4:6,:);     % 1/m

z = -H/2:10*1e-6:H/2;                     % z for whole laminate
for i = 1:nots
    eps(:,:,i) = repmat(eps0(:,i), 1, length(z)) + z .* kappa(:,i);
end
%% Stresses
sigma = zeros(size(eps));
sigma_loc = zeros(size(eps));
for i = 1:length(z)

    for j = 1:length(h)-1
        if z(i)>=h(j) && z(i)<h(j+1)
             ply = j;
             break
        end
    end
    if z(i)==h(end)
        ply = length(h)-1;
    end

    for j = 1:nots
        sigma(:,i,j) = Qbar(:,:,ply) * eps(:,i,j);
    end
    % Local stresses
    s = sin(theta(ply)*pi/180);
    c = cos(theta(ply)*pi/180);
    % Transformation matrix
    T = [c^2    s^2     2*s*c;
         s^2    c^2     -2*s*c;
         -s*c   s*c     c^2-s^2];  % -
    for j = 1:nots
        sigma_loc(:,i,j) = T * sigma(:,i,j);
    end
    
end

%% Strength Ratio

SR = zeros(length(z),nots);
% Tsai-Wu Criterion
H1 = 1/sigma_1_T_ult - 1/sigma_1_C_ult;     % 1/Pa
H11 = 1/(sigma_1_T_ult*sigma_1_C_ult);      % 1/Pa^2
H2 = 1/sigma_2_T_ult - 1/sigma_2_C_ult;     % 
H22 = 1/(sigma_2_T_ult*sigma_2_C_ult);
H6 = 0;
H66 = 1/tau_12_ult^2;
H12 = -1/2 * sqrt(1/(sigma_1_T_ult*sigma_1_C_ult*sigma_2_T_ult*sigma_2_C_ult)); % Mises-Hencky Criterion
for i=1:length(z)
    for j = 1:nots
        p = [H11*sigma_loc(1,i,j)^2+H22*sigma_loc(2,i,j)^2+H66*sigma_loc(3,i,j)^2+...
             H12*sigma_loc(1,i,j)*sigma_loc(2,i,j)...
             H1*sigma_loc(1,i,j)+H2*sigma_loc(2,i,j)+H6*sigma_loc(3,i,j) -1];
        SR(i,j) = max(roots(p));
    end
end

SR_tw = SR;

%Tsai-Hill Criterion
SR = zeros(length(z),nots);
for i=1:length(z)
    for j = 1:nots
        if sigma_loc(1,i,j) < 0 X1 = sigma_1_C_ult; else X1 = sigma_1_T_ult; end
        if sigma_loc(2,i,j) < 0 X2 = sigma_2_C_ult; else X2 = sigma_2_T_ult; end
        FI = sigma_loc(1,i,j)^2/X1^2 - sigma_loc(1,i,j)*sigma_loc(2,i,j)/X1^2 + sigma_loc(2,i,j)^2/X2^2 + sigma_loc(3,i,j)^2/tau_12_ult^2;
        SR(i,j) = 1/sqrt(FI);
    end
end

SR_th = SR;

%Max Stress Criterion
FI_1 = reshape(max(sigma_loc(1,:,:),0)/sigma_1_T_ult,[length(z),nots]);
FI_2 = reshape(max(-sigma_loc(1,:,:),0)/sigma_1_C_ult,[length(z),nots]);
FI_3 = reshape(max(sigma_loc(2,:,:),0)/sigma_2_T_ult,[length(z),nots]);
FI_4 = reshape(max(-sigma_loc(2,:,:),0)/sigma_2_C_ult,[length(z),nots]);
FI_5 = reshape(abs(sigma_loc(3,:,:))/tau_12_ult,[length(z),nots]);

SR_ms = 1./max(cat(3,FI_1, FI_2, FI_3, FI_4, FI_5),[],3);
%disp([SR_ms, SR_th, SR_tw])
SR = min(cat(3,SR_ms, SR_th, SR_tw),[],3);
%% Most critical case
[M1,I1]=min(SR,[],1);
[M2,I2]=min(M1);
i = I1(I2);
j = I2;
min_SR = SR(i,j);
s1 = sigma_loc(1,i,j);
s2 = sigma_loc(2,i,j);
t12 = sigma_loc(3,i,j);
%% Output
fprintf('Total mass: %.2f g\n',mass*1e3)
fprintf('Minimum strength ratio:\t\t\t\t\t%.2f\n',min(SR,[],'all'))
fprintf('Minimum strength ratio (Max Stress):\t%.2f\n',min(SR_ms,[],'all'))
fprintf('Minimum strength ratio (Tsai-Wu):\t\t%.2f\n',min(SR_tw,[],'all'))
fprintf('Minimum strength ratio (Tsai-Hill):\t\t%.2f\n',min(SR_th,[],'all'))
toc
%% Plotting Strain and Stress
step = j;

f1 = figure('name','Strain','numberTitle','off');
f1.Position = [403  20   560   670];
grid on;
plot(eps(:,:,step),z*1e3,LineWidth=1.5)
grid on;
xlabel('Strain');
ylabel('z (mm)');
legend('\epsilon_x','\epsilon_y','\gamma_{xy}',Location='best')
set(gca, 'YDir','reverse')
yticks(h*1e3)
y_tick_labels = num2cell(h);
for i=1:length(y_tick_labels)
    if mod(i,2) == 0
        y_tick_labels(i)={[]};
    else
        y_tick_labels(i)={num2str(h(i)*1e3,'%.2f')};
    end
end
yticklabels(y_tick_labels);
% set(gca, 'FontSize', 10)

exportgraphics(f1,'Plots/strain_laminate.eps', BackgroundColor='none',ContentType='vector')
%%
f2 = figure('name','Stress','numberTitle','off');
f2.Position = [403  20   560   670];
plot(sigma(:,:,step)*1e-6,z*1e3,LineWidth=1.5)
grid on;
legend('\sigma_x','\sigma_y','\tau_{xy}',Location='best')
set(gca, 'YDir','reverse')
yticks(h*1e3)
y_tick_labels = num2cell(h);
for i=1:length(y_tick_labels)
    if mod(i,2) == 0
        y_tick_labels(i)={[]};
    else
        y_tick_labels(i)={num2str(h(i)*1e3,'%.2f')};
    end
end
yticklabels(y_tick_labels);
ylabel('z (mm)')
xlabel('\sigma (MPa)')
exportgraphics(f2,'Plots/stress_laminate.eps', BackgroundColor='none',ContentType='vector')
%%
f3 = figure('name','Strength Ratio','numberTitle','off');
f3.Position = [403  20   560   670];
plot(SR(:,step),z*1e3,LineWidth=1.5)
grid on;
set(gca, 'YDir','reverse')
set(gca, 'XScale', 'log')
yticks(h*1e3)
y_tick_labels = num2cell(h);
for i=1:length(y_tick_labels)
    if mod(i,2) == 0
        y_tick_labels(i)={[]};
    else
        y_tick_labels(i)={num2str(h(i)*1e3,'%.2f')};
    end
end
yticklabels(y_tick_labels);
ylabel('z (mm)')
xlabel('Strength Ratio')
exportgraphics(f3,'Plots/sr_laminate.eps', BackgroundColor='none',ContentType='vector')
%%
f4 = figure('name','Strength Ratio','numberTitle','off');
plot(100*(1:nots)/nots,min(SR),LineWidth=1.5);
xlabel('Percentage of stance gait (%)')
grid on;
% set(gca, 'YScale', 'log')
ylabel('Strength Ratio')
exportgraphics(f4,'Plots/sr_time.eps', BackgroundColor='none',ContentType='vector')
%% Failure Loci sigma1 sigma2
f5 = figure('name','Failure Loci sigma1 sigma2','numberTitle','off');
hold on;
grid on;
xlabel('\sigma_1 (MPa)');
ylabel('\sigma_2 (MPa)');
plot([-sigma_1_C_ult sigma_1_T_ult sigma_1_T_ult -sigma_1_C_ult -sigma_1_C_ult], ...
     [-sigma_2_C_ult -sigma_2_C_ult sigma_2_T_ult sigma_2_T_ult -sigma_2_C_ult],...
     DisplayName='Max Stress',LineWidth=1.5,Color='#1b9e77');
% Tsai Wu
f_tw = @(s1,s2) H1.*s1 + H2.*s2 + H11.*s1.^2 + H22.*s2.^2 + H12.*s1.*s2 +H6.*t12 + H66.*t12.^2 - 1;
fimplicit(f_tw,[-sigma_1_C_ult sigma_1_T_ult -sigma_2_C_ult sigma_2_T_ult]*2,...
    DisplayName='Tsai-Wu',LineWidth=1.5,Color='#d95f02');
% Tsai Hill
f_th_11 = @(s1,s2) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + ...
    s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
f_th_12 = @(s1,s2) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + ...
    s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
f_th_21 = @(s1,s2) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + ...
    s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
f_th_22 = @(s1,s2) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + ...
    s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
fimplicit(f_th_11,[0 sigma_1_T_ult 0 sigma_2_T_ult]*2,'r',...
    DisplayName='Tsai-Hill',LineWidth=1.5,Color='#7570b3');
fimplicit(f_th_12,[0 sigma_1_T_ult -sigma_2_C_ult 0]*2,'r',DisplayName='',...
    LineWidth=1.5,HandleVisibility='off',Color='#7570b3');
fimplicit(f_th_21,[-sigma_1_C_ult 0 0 sigma_2_T_ult]*2,'r',DisplayName='',...
    LineWidth=1.5,HandleVisibility='off',Color='#7570b3');
fimplicit(f_th_22,[-sigma_1_C_ult 0 -sigma_2_C_ult 0]*2,'r',DisplayName='',...
    LineWidth=1.5,HandleVisibility='off',Color='#7570b3');
% Our Design
scatter(s1,s2,"filled",DisplayName='Our Design',MarkerFaceColor='#a2142f');
title(['\tau_{12} = ' num2str(round(t12*1e-6,2)) ' MPa']);
xlim([-3000 3000]*1e6);
ylim([-250 150]*1e6);
h=gca;
h.XTickLabel = h.XTick * 1e-6;
h.YTickLabel = h.YTick * 1e-6; 
% axis equal;
legend;
exportgraphics(f5,'Plots/loci_s1_s2.eps', BackgroundColor='none',ContentType='vector')
%% sigma1 tau12
f6 = figure('name','Failure Loci sigma1 tau12','numberTitle','off');
hold on;
grid on;
xlabel('\sigma_1 (MPa)');
ylabel('\tau_{12} (MPa)');
plot([-sigma_1_C_ult sigma_1_T_ult sigma_1_T_ult -sigma_1_C_ult -sigma_1_C_ult], ...
     [-tau_12_ult -tau_12_ult tau_12_ult tau_12_ult -tau_12_ult],...
     DisplayName='Max Stress',LineWidth=1.5,Color='#1b9e77');
% Tsai Wu
f_tw = @(s1,t12) H1.*s1 + H2.*s2 + H11.*s1.^2 + H22.*s2.^2 + H12.*s1.*s2 +H6.*t12 + H66.*t12.^2 - 1;
fimplicit(f_tw,[-sigma_1_C_ult sigma_1_T_ult -tau_12_ult tau_12_ult]*2,...
    DisplayName='Tsai-Wu',LineWidth=1.5,Color='#d95f02');
% Tsai Hill
if s2 > 0
    f_th_1 = @(s1,t12) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + ...
        s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
    f_th_2 = @(s1,t12) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + ...
        s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
else
    f_th_1 = @(s1,t12) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + ...
        s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
    f_th_2 = @(s1,t12) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + ...
        s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
end
fimplicit(f_th_1,[0 sigma_1_T_ult -tau_12_ult tau_12_ult]*2,'r',...
    DisplayName='Tsai-Hill',LineWidth=1.5,Color='#7570b3');
fimplicit(f_th_2,[-sigma_1_C_ult 0 -tau_12_ult tau_12_ult]*2,'r',...
    DisplayName='',LineWidth=1.5,HandleVisibility='off',Color='#7570b3');
% Our Design
scatter(s1,t12,"filled",DisplayName='Our Design',MarkerFaceColor='#a2142f');
title(['\sigma_{2} = ' num2str(round(s2*1e-6,2)) ' MPa']);
xlim([-2500 3500]*1e6);
% ylim([-250 150]*1e6);
h=gca;
h.XTickLabel = h.XTick * 1e-6;
h.YTickLabel = h.YTick * 1e-6; 
% axis equal;
legend
exportgraphics(f6,'Plots/loci_s1_t12.eps', BackgroundColor='none',ContentType='vector')
%% sigma2 tau12
f7 = figure('name','Failure Loci sigma2 tau12','numberTitle','off');
hold on;
grid on;
xlabel('\sigma_2 (MPa)');
ylabel('\tau_{12} (MPa)');
plot([-sigma_2_C_ult sigma_2_T_ult sigma_2_T_ult -sigma_2_C_ult -sigma_2_C_ult], ...
     [-tau_12_ult -tau_12_ult tau_12_ult tau_12_ult -tau_12_ult],...
     DisplayName='Max Stress',LineWidth=1.5,Color='#1b9e77');
% Tsai Wu
f_tw = @(s2,t12) H1.*s1 + H2.*s2 + H11.*s1.^2 + H22.*s2.^2 + H12.*s1.*s2 +H6.*t12 + H66.*t12.^2 - 1;
fimplicit(f_tw,[-sigma_2_C_ult sigma_2_T_ult -tau_12_ult tau_12_ult]*2,...
    DisplayName='Tsai-Wu',LineWidth=1.5,Color='#d95f02');
% Tsai Hill
if s1 > 0
    f_th_1 = @(s2,t12) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
    f_th_2 = @(s2,t12) s1.^2/sigma_1_T_ult^2 - s1.*s2./sigma_1_T_ult^2 + s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
else
    f_th_1 = @(s2,t12) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + s2.^2./sigma_2_T_ult^2 + t12.^2./tau_12_ult^2 - 1;
    f_th_2 = @(s2,t12) s1.^2/sigma_1_C_ult^2 - s1.*s2./sigma_1_C_ult^2 + s2.^2./sigma_2_C_ult^2 + t12.^2./tau_12_ult^2 - 1;
end
fimplicit(f_th_1,[0 sigma_2_T_ult -tau_12_ult tau_12_ult]*2,'r',...
    DisplayName='Tsai-Hill',LineWidth=1.5,Color='#7570b3');
fimplicit(f_th_2,[-sigma_2_C_ult 0 -tau_12_ult tau_12_ult]*2,'r',...
    DisplayName='',LineWidth=1.5,HandleVisibility='off',Color='#7570b3');
% Our design
scatter(s2,t12,"filled",DisplayName='Our Design',MarkerFaceColor='#a2142f');
title(['\sigma_{1} = ' num2str(round(s1*1e-6,2)) ' MPa']);
xlim([-250 100]*1e6);
h=gca;
h.XTickLabel = h.XTick * 1e-6;
h.YTickLabel = h.YTick * 1e-6; 
% axis equal;
legend;
exportgraphics(f7,'Plots/loci_s2_t12.eps', BackgroundColor='none',ContentType='vector')