%% Data Indexing

close all;
clearvars;
clc;
data_files = dir('output_*y.mat');

%% Full Development

% Highest refinement selection
load(data_files(end).name);
delta = max(Y_N);

% Variable profiles vs X-axis
figure(1);
plot((X_E(3:end) - 0.1) / delta, V1(3:end, ceil(NY/2)));
xlabel("x (deltas)");
ylabel("V [m/s]");

figure(2);
DUDX = filter([1/2, 0, -1/2], 1, U1, [], 1) ./ filter([1/2, 0, -1/2], 1, repmat(X_E', 1, NY));
plot((X_E(3:end-1) - 0.1) / delta, DUDX(3:end-1, ceil(NY/2)));
xlabel("x (deltas)");
ylabel("DUDX [1/s]");

%% Grid Independence

for i=1:length(data_files)
    % Import Data
    load(data_files(i).name);
    x_study = ceil(0.8 * NX);
    y_study = ceil(NY/2);
    delta = max(Y_N);
    
    % X-Velocity Y-Profile
    figure(3);
    plot(Y_C, U1(x_study, :));
    hold on;

    % Y-Velocity Y-Profile
    figure(4);
    plot(Y_N, V1(x_study, :));
    hold on;

    % Pressure X-Profile
    figure(5);
    plot((X_C(3:end-1) - 0.1) / delta, P1(3:end-1, y_study));
    hold on;

    fprintf("- %i x %i: YPLUS > %d\n", NX, NY, min(YPLS(3:end,1)));
end

figure(3);
xlabel('y [m]');
ylabel('u [m/s]');
legend('250 x 5', '500 x 10', '750 x 15', '1000 x 20', '1500 x 30');
hold off;

figure(4);
xlabel('y [m]');
ylabel('v [m/s]');
legend('250 x 5', '500 x 10', '750 x 15', '1000 x 20', '1500 x 30');
hold off;

figure(5);
xlabel('x [deltas]');
ylabel('p [Pa]');
legend('250 x 5', '500 x 10', '750 x 15', '1000 x 20', '1500 x 30');
hold off;

%% Stress Study

load(data_files(end).name);
rho = 9.9823E2;
nu = 1E-6;
mu = rho * nu;
x_study = ceil(0.8 * NX);

tau_wall= STRS(x_study, 1) * rho;
tau_vis= mu * DUDY(x_study, :); 
tau = (1 - (Y_C / delta)) * tau_wall;
tau_re= tau - tau_vis;

figure(6);
plot(Y_C / delta, tau_re / tau_wall);
hold on;
plot(Y_C / delta, tau_vis / tau_wall);
plot(Y_C / delta, tau / tau_wall);
xlabel('y [deltas]');
ylabel('Stress [normalized]');
legend('TauRe / TauWall', 'TauVis / TauWall', 'Tau / TauWall');
hold off;

%% Axial Velocity Profile

U0 = U1(end-1, end);

figure(7);
plot(Y_C / delta, U1(x_study, :) / U0);
hold on;
plot(Y_C / delta, Y_C / delta .* (2 - Y_C / delta))
hold off;
xlabel('y/delta [normalized]');
ylabel('u/U0 [normalized]');
legend(["Simulated Turbulent Profile"; "Analytic Laminar Profile"],...
        'Location','best');

%% Axial Velocity vs Model

% Simulated
u_tau = sqrt(tau_wall / rho);
Yp = Y_C * u_tau / nu;
Up = U1(x_study,:) / u_tau;

% Calculated
% y_pls_calc = logspace(log10(30), log10(u_tau * (0.1 * delta / nu)), 1e2);
y_pls_calc = Yp;
k = 0.41;
E = 8.6;
u_plus_calc = log(E * y_pls_calc) / k;

% Plot
figure(8);
semilogx(Yp, Up, 'marker', '*');
hold on;
semilogx(y_pls_calc, u_plus_calc, 'marker', '*');
xlabel("y+");
ylabel("u+");
semilogx([u_tau * (0.1 * delta / nu), u_tau * (0.1 * delta / nu)], [min(Up), max(Up)]); % [30, y+(y=delta/10)]
legend('Simulated', 'Calculated', 'Model limit of operation');
hold off;

%% Skin Friction Coefficient

% Parameters
c_f_sim = 2 * tau_wall / (rho * U0 * U0);
UB = 1;
Re = 2 * delta * UB / nu;

% Expressions of the relation
c_f_range = logspace(log10(c_f_sim / 10), log10(c_f_sim * 10), 1e2);
c_f_left = sqrt(2 ./ c_f_range);
c_f_right = log(Re ./ (sqrt(8 ./ c_f_range) - 4.8)) / k + 5.2 + 0.5;

% Plot
figure(9);
semilogx(c_f_range, c_f_left);
hold on;
semilogx(c_f_range, c_f_right);
scatter(c_f_sim,...
    (sqrt(2 ./ c_f_sim) + log(Re ./ (sqrt(8 ./ c_f_sim) - 4.8)) / k + 5.2 + 0.5) / 2,...
    50, 'LineWidth', 3);
xlim([c_f_sim/10, c_f_sim * 10]);
xlabel("C_f");
legend("LHS", "RHS", "Calculated Value");
hold off;

%% Reynolds variation

data_files = dir('output_Re*.mat');
Res = [1e4 2e4 5e4 8e4 1e5];
cfs_sim = zeros(size(data_files));
cfs_calc = zeros(size(data_files));

rho = 9.9823E2;
syms cf;
for i=1:length(data_files)
    % Import Data
    load(data_files(i).name);
    U0 = U1(end-1, end);
    x_study = ceil(0.8 * NX);

    tau_wall= STRS(x_study, 1) * rho;
    cfs_sim(i) = 2 * tau_wall / (rho * U0 * U0);

    eqn = sqrt(2/cf) == log(Res(i)*(sqrt(8/cf)-4.8)^-1)/0.41+5.2+0.5;
    cfs_calc(i) = vpasolve(eqn,cf,cfs_sim(i));
end

figure(10);
semilogx(Res, cfs_sim, 'Marker','*');
hold on;
semilogx(Res, cfs_calc, 'Marker','*');
xlabel("Re");
ylabel("C_f");
hold off;

%% Cf vs Ypls

data_files = dir('output_yp*.mat');
yps = zeros(size(data_files));
cfs = zeros(size(data_files));

rho = 9.9823E2;
for i=1:length(data_files)
    load(data_files(i).name);
    x_study = ceil(0.8 * NX);
    U0 = U1(end-1, end);

    tau_wall= STRS(x_study, 1) * rho;
    yps(i) = min(YPLS(3:end,1));
    cfs(i) = 2 * tau_wall / (rho * U0 * U0);
end

figure(11);
semilogx(yps,cfs,'Marker','*');
xlabel("Y_+");
ylabel("C_f");