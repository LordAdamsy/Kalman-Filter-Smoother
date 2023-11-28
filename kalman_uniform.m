clc;clear;close all;

load('vel_data.mat');
load('pos_data.mat');

% parameters

L = length(pos_e);
z_pos = [pos_e; pos_n];
z_pos_vel = [pos_e; pos_n; vel_e; vel_n];
t = 5e-3;
F = [1, 0, t, 0; 0, 1, 0, t; 0, 0, 1, 0; 0, 0, 0, 1];
C_pos = [1, 0, 0, 0; 0, 1, 0, 0];
C_pos_vel = eye(4);
Rn_pos = eye(2)*36;
Rn_pos_vel = [eye(2)*36, zeros(2, 2); zeros(2, 2), eye(2)];
Q = [t^4/4, 0, t^3/2, 0; 0, t^4/4, 0, t^3/2; t^3/2, 0, t^2, 0; 0, t^3/2, 0, t^2];

s_pos = zeros(4, L);
s_pos_vel = zeros(4, L);
x_pos = ones(4, 1)*10;
x_pos_vel = ones(4, 1)*10;
K_pos = eye(4)*1000;
K_pos_vel = eye(4)*1000;

% Kalman filter
for i = 1 : L

    % with postion infomation
    x_n_n_1_pos = F * x_pos;
    alpha_pos = z_pos(:, i) - C_pos * x_n_n_1_pos;
    K_n_n_1_pos = F * K_pos * F' + Q;
    R_pos = C_pos * K_n_n_1_pos * C_pos' + Rn_pos;
    Gf_pos = K_n_n_1_pos * C_pos' / R_pos;
    x_pos = x_n_n_1_pos + Gf_pos * alpha_pos;
    K_pos = (eye(4) - Gf_pos * C_pos) * K_n_n_1_pos;
    s_pos(:, i) = x_pos;

    % with position and velocity infomation
    x_n_n_1_pos_vel = F * x_pos_vel;
    alpha_pos_vel = z_pos_vel(:, i) - C_pos_vel * x_n_n_1_pos_vel;
    K_n_n_1_pos_vel = F * K_pos_vel * F' + Q;
    R_pos_vel = C_pos_vel * K_n_n_1_pos_vel * C_pos_vel' + Rn_pos_vel;
    Gf_pos_vel = K_n_n_1_pos_vel * C_pos_vel' / R_pos_vel;
    x_pos_vel = x_n_n_1_pos_vel + Gf_pos_vel * alpha_pos_vel;
    K_pos_vel = (eye(4) - Gf_pos_vel * C_pos_vel) * K_n_n_1_pos_vel;
    s_pos_vel(:, i) = x_pos_vel;
end

figure
plot(s_pos(1, :), s_pos(2, :), 'LineWidth', 1.3);
hold on;
plot(s_pos_vel(1, :), s_pos_vel(2, :), 'LineWidth', 1.3);
hold on;
plot(pos_e, pos_n, 'LineWidth', 1.3);
hold on;
xlabel('east position');
ylabel('north position');
legend('pos', 'pos and vel', 'observed');
title('kalman filter with uniform model');
grid on;