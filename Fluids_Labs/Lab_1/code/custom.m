%% Full-Development Analysis

clc;
close all;
load('output.mat');
figure();
for i=6:2:IX-2
    plot(Y_N, U1(i, :), 'DisplayName', string(X_E(i)));
    hold on;
end
legend(strcat(string(2 * (X_E(6:2:end-2) - 0.01) / Y_N(end)), ' deltas'));
% title('X-velocity Y-profile per delta-step');
xlabel('y (m)');
ylabel('u (m/s)');
hold off;

%% Full-Development Analysis (3D Surface)

% delta = 0.005;
% figure();
% surf(repmat(X_E' ./ delta, [1, IY]), repmat(Y_N ./ delta, [IX, 1]), U1);
% title('X-Velocity Y-profile');
% xlabel('x (normalized by delta)');
% ylabel('y (normalized by delta)');
% zlabel('u (m/s)');

%% Grid-Independence Analysis

% Data extracted with Phoenics plot generator

b_coeffs = [-1/280, 4/105, -1/5, 4/5, 0, -4/5, 1/5, -4/105, 1/280];
% b_coeefs = [1/2, 0, -1/2];
% b_coeffs = [1, -1];

figure();

csv_list_u_x = dir('Results/velocity_ind_x/X-Velocity_profile_*.csv');
csv_list_p_x = dir('Results/pressure_ind_x/Pressure_profile_*.csv');
n_csv_x = length(csv_list_u_x);
u_x_data = cell(1, n_csv_x);
p_x_data = cell(1, n_csv_x);
x_analysis_n = 10:5:50;
bulk_u_x = zeros(size(u_x_data));
p_grad_x = zeros(size(p_x_data));

subplot(2,1,1);
u_x_legend = strings(size(x_analysis_n));
for i=1:n_csv_x
    u_x_data{i} = readmatrix(strcat(csv_list_u_x(i).folder, '/', csv_list_u_x(i).name));
    p_x_data{i} = readmatrix(strcat(csv_list_p_x(i).folder, '/', csv_list_p_x(i).name));
    bulk_u_x(i) = sum(u_x_data{i}(:,2)...
        ./ length(u_x_data{i}));

    u_x_data{i} = [u_x_data{i},...
        filter(b_coeffs, 1, u_x_data{i}(:,2)) ./ filter([1 -1], 1, u_x_data{i}(:,1))];
    u_x_legend(i) = strcat(num2str(x_analysis_n(i)), " cells");
    plot(u_x_data{i}(length(b_coeffs) + 1:end,1),...
        u_x_data{i}(length(b_coeffs) + 1:end,3));
    hold on;

    p_grad_x(i) = median(filter(b_coeffs, 1, p_x_data{i}(:,2))...
        ./ filter([1 -1], 1, p_x_data{i}(:,1)));
end
legend(u_x_legend);
title('X-velocity derivative Y-profile per X-refinement');
xlabel('y (m)');
ylabel('du/dy (1/s)');
xlim([0 Y_N(end)]);
hold off;

csv_list_u_y = dir('Results/velocity_ind_y/X-Velocity_profile_*.csv');
csv_list_p_y = dir('Results/pressure_ind_y/Pressure_profile_*.csv');
n_csv_y = length(csv_list_u_y);
u_y_data = cell(1, n_csv_y);
p_y_data = cell(1, n_csv_y);
y_analysis_n = 5:5:50;
bulk_u_y = zeros(size(u_y_data));
p_grad_y = zeros(size(p_y_data));

subplot(2,1,2);
u_y_legend = strings(size(y_analysis_n));
for i=1:n_csv_y
    u_y_data{i} = readmatrix(strcat(csv_list_u_y(i).folder, '/', csv_list_u_y(i).name));
    p_y_data{i} = readmatrix(strcat(csv_list_p_y(i).folder, '/', csv_list_p_y(i).name));
    bulk_u_y(i) = sum(u_y_data{i}(:,2)...
        ./ length(u_y_data{i}));

    u_y_data{i} = [u_y_data{i},...
        filter(b_coeffs, 1, u_y_data{i}(:,2)) ./ filter([1 -1], 1, u_y_data{i}(:,1))];
    u_y_legend(i) = strcat(num2str(y_analysis_n(i)), " cells");
    plot(u_y_data{i}(length(b_coeffs):end,1),...
        u_y_data{i}(length(b_coeffs):end,3));
    hold on;

    p_grad_y(i) = median(filter(b_coeffs, 1, p_y_data{i}(:,2))...
        ./ filter([1 -1], 1, p_y_data{i}(:,1)));
end
legend(u_y_legend);
title('X-velocity derivative Y-profile per Y-refinement');
xlabel('y (m)');
ylabel('du/dx (1/s)');
xlim([0 Y_N(end)]);
hold off;

figure();
subplot(2,1,1);
plot(x_analysis_n, p_grad_x);
title('P-gradient vs X-refinement');
subplot(2,1,2);
plot(y_analysis_n, p_grad_y);
title('P-gradient vs Y-refinement');

%% Analytic Validation (using 40x40)

delta = 5e-3;
mu = 1.006e-3;
Ub = 5e-3;

dpedx_calc = - 3 * mu * Ub / (delta^2);
dpedx_sim = - sqrt(p_grad_x(7) * p_grad_y(8));
dpedx_err = abs(dpedx_sim - dpedx_calc) / abs(dpedx_calc);
fprintf(['The expected value for the X-derivative of the pressure was '...
    '%d Pa/m,\nwhereas the simulated value was %d Pa/m,\nthus leading to a '...
    'relative error of %d.\n'], dpedx_calc, dpedx_sim, dpedx_err);

% The first and last 10 samples are removed 
n_removed = 50;
y = u_x_data{7}(n_removed:end-n_removed,1);
u_calc = - (delta / (2 * mu)) * dpedx_calc * y .* (2 - y / delta);
u_sim = sqrt(u_x_data{7}(n_removed:end-n_removed, 2) .* u_y_data{8}(n_removed:end-n_removed, 2));
u_calc(u_calc==0) = u_sim(u_calc==0);
[u_err_Linf, max_ind] = max(abs((u_calc - u_sim) ./ u_calc));
u_err_L2 = sqrt(sum(((u_calc - u_sim) ./ u_calc).^2) / length(y));
figure();
plot(y, u_calc);
hold on;
plot(y, u_sim);
title("X-velocity profiles on the fully-developed region");
legend(["Analytical", "Simulated"]);
xlabel("y (m)");
ylabel("u (m/s)");
xlim([y(1), y(end)]);
hold off;
fprintf(['The maximum relative error between the simulated X-velocity and ' ...
    'the analytical one was of %d, found at y = %d m.\n'], u_err_Linf, y(max_ind));
fprintf(['The average root-mean-squared relative error between the simulated X-velocity and ' ...
    'the analytical one was of %d.\n'], u_err_L2);

tau_calc = dpedx_calc * (delta - y);
tau_sim = - mu * 0.5 *(u_x_data{7}(n_removed:end-n_removed, 3) + u_y_data{8}(n_removed:end-n_removed, 3));
figure();
plot(y, tau_calc);
hold on;
plot(y, tau_sim);
title("Shear-stress profiles on the fully-developed region");
legend(["Analytical", "Simulated"]);
xlabel("y (m)");
ylabel("τ (Pa)");
xlim([y(1), y(end)]);
hold off;

%% Vorticity

% b_coeffs = [-1/12, 2/3, 0, -2/3, 1/12];
% DVDX = filter(b_coeffs, 1, V1, [], 1) ./ repmat(filter([1 -1], 1, X_C)', 1, IY);
% DUDY = - filter(b_coeffs, 1, U1, [], 2) ./ repmat(filter([1 -1], 1, Y_C), IX, 1);
% Vort = abs(DVDX - DUDY);
% 
% figure();
% surf(repmat((X_C(5:end-4))', [1, IY-8]), repmat(Y_C(5:end-4), [IX-8, 1]), VOR1(5:end-4, 5:end-4));
% title('Vorticity profile');
% xlabel('x (m)');
% ylabel('y (m)');
% zlabel('Vorticity (s⁻¹)');

%% Vorticity (Francesco's data)

% clearvars;
delta = 5e-3;
load("Results/k_x.mat");
load("Results/k_y.mat");
load("Results/vort.mat");
figure();
surf(repmat(k_x', [1, length(k_y)]), repmat(k_y, [length(k_x), 1]), VOR1);
title('Vorticity profile');
xlabel('x (m)');
xlim([k_x(1) k_x(end)]);
ylabel('y (m)');
ylim([0 delta]);
zlabel('Vorticity (s⁻¹)');