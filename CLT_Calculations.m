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
% [45/-45/45/-45/0/90/0/90/-45/0/45/90]_s
% theta = [45 -45 45 -45 0 90 0 90 -45 0 45 90];
theta = 2*ones(1,8)*45;    % degree
if core == true
    theta = [theta 0 flip(theta)];                       % degree (Symmetric)
    n = size(theta,2);  % number of plies
    H = n*t+t_core;        % m % Total width of the lamimate
    for i = 0:length(theta)/2
        h(i+1) = -H/2 + i*t; % m
    end
    h = [h -flip(h)];
else
    theta = [theta flip(theta)];                       % degree (Symmetric)
    n = size(theta,2);  % number of plies
    H = n*t;        % m % Total width of the lamimate
    for i = 0:length(theta)
        h(i+1) = -H/2 + i*t; % m
    end
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
shoe_size = 42;     % eu
L_data = 230e-3;    % m % Total Lenght of the foot
% L_model = 206e-3;  % m for proted design
L_model = ((shoe_size - 2 ) * 20 / 3)*1e-3;
b_rear = 0.225 * L_model;  % Lenght of the front part of the foot
b_front = L_model - b_rear;    % Lenght of the rear part of the foot
a = 0.31 * L_model; % Width of the foot
delta = 7.5;    % deg % the angle between foot axis and the walking direction
Iy = (1/12) * a*H^3; % m^4
Iz = (1/12) * H*a^3; % m^4
J = Iy + Iz; % m^4

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
design_mass = 130;

number_of_time_steps = length(spi);
nots = number_of_time_steps;
% Toe contact
Fz = -F_foot_ground_yp(spi)*design_mass/mass_data;
Fy = F_foot_ground_xp(spi) * sind(delta)*design_mass/mass_data;
Fx = F_foot_ground_xp(spi) * cosd(delta)*design_mass/mass_data;
b = CoP_xp(spi)*1e-3 * L_model / L_data - b_rear;
% b = abs(b);
%% Loadings
    
Nx = Fx/(H*a)*H + Fy.*b*(a/2)*H/Iz;    % N/m
Ny = Fy./(H*b)*H;    % N/m
Nxy = Fy/(H*a)*H;  % N/m
N = [Nx Ny Nxy]';

Mx = zeros(size(b));         % N*m/m
My = Fz.*b*(H^3/12)/Iy;      % N*m/m
Mxy = -Fz*a^4/(16*J)*(-1/4*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2)))) + 1/2*H/a*(H^2/a^2 + 1)^(3/2));                 % N*m/m
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
% Mises-Hencky Criterion
H12 = -1/2 * sqrt(1/(sigma_1_T_ult*sigma_1_C_ult*sigma_2_T_ult*sigma_2_C_ult));
for i=1:length(z)
    for j = 1:nots
        p = [H11*sigma_loc(1,i,j)^2+H22*sigma_loc(2,i,j)^2+H66*sigma_loc(3,i,j)^2+...
             H12*sigma_loc(1,i,j)*sigma_loc(2,i,j)...
             H1*sigma_loc(1,i,j)+H2*sigma_loc(2,i,j)+H6*sigma_loc(3,i,j) -1];
        SR(i,j) = max(roots(p));
    end
end
%% Output
fprintf('Total mass: %.2f g\n',mass*1e3)
fprintf('Minimum strength ratio: %.2f\n',min(min(SR)))
toc
%% Plotting Strain and Stress
f1 = figure('name','Strain','numberTitle','off');
grid on;

plot(eps(:,:,35),z*1e3,LineWidth=1.5)
grid on;
legend('\epsilon_x','\epsilon_y','\epsilon_{xy}',Location='best')
set(gca, 'YDir','reverse')
yticks(h*1e3)



figure('name','Stress','numberTitle','off');
plot(sigma(:,:,35)*1e-6,z*1e3,LineWidth=1.5)
grid on;
legend('\sigma_x','\sigma_y','\tau_{xy}',Location='best')
set(gca, 'YDir','reverse')
yticks(h*1e3)
ylabel('z (mm)')
xlabel('\sigma (MPa)')

figure('name','Strength Ratio','numberTitle','off');
plot(SR(:,35),z*1e3,LineWidth=1.5)
grid on;
set(gca, 'YDir','reverse')
set(gca, 'XScale', 'log')
yticks(h*1e3)
ylabel('z (mm)')
xlabel('Strength Ratio')

figure('name','Strength Ratio','numberTitle','off');
plot(min(SR),LineWidth=1.5)
grid on;
set(gca, 'YScale', 'log')
ylabel('Strength Ratio')