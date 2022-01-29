%% Data Indexing

close all;
clearvars;
clc;

files = dir("*.mat");

%% Grid Independence

titles = ["30x15"; "50x25"; "80x40"; "100x50"];
for i=1:length(files)
    load(files(i).name);

    figure(1);
    plot(Y_C, U1(end-1,:));
    hold on;

    figure(2);
    plot(Y_N, V1(end-1, :));
    hold on;
end

figure(1);
xlabel("y [m]");
ylabel("U [m/s]");
legend(titles);
hold off;

figure(2);
xlabel("y [m]");
ylabel("V [m/s]");
hold off;
legend(titles);

%% Blasius Velocity

U_inf = 1;
nu = 1.5445e-5;

lim_inf = NX/2;
lim_sup = NX - 1;

titles = cell(1, length(lim_inf:15:lim_sup) + 1);
index = 1;
for j = lim_inf:15:lim_sup
    titles{index} = num2str(X_E(j) - 0.1);
    eta = Y_C * sqrt(U_inf ./ (nu * (X_E(j) - 0.1)));
    U_norm = U1(j, :) / U_inf;
    V_norm = sqrt(U_inf * (X_E(j) - 0.1) / nu) * V1(j, :) / U_inf;

    figure(3);
    plot(eta(eta<=6), U_norm(eta<=6));
    hold on;

    figure(4);
    plot(eta(eta<=6), V_norm(eta<=6));
    hold on;

    index = index + 1;
end
titles{end} = "reference";

eta_ref = [0, 0.4, 1, 1.4, 2, 2.4, 3, 3.4, 4, 4.4, 5, 6, 7];
f_ref = [0, 0.026, 0.1655, 0.323, 0.65, 0.992, 1.397,...
    1.7469, 2.305, 2.692, 3.283, 4.279, 5.279];
f_p_ref = [0, 0.133, 0.3298, 0.4563, 0.63, 0.7289, 0.846,...
    0.9017, 0.956, 0.9758, 0.992, 0.998, 0.999];

eta = eta(eta<=6);

figure(3);
plot(eta_ref, f_p_ref, '--');
xlabel("eta");
ylabel("U/U_inf");
xlim([eta(1), eta(end)]);
legend(titles, 'location', 'best');
hold off;

figure(4);
plot(eta_ref, (eta_ref .* f_p_ref - f_ref) / 2, '--');
xlabel("eta");
ylabel("sqrt(U_inf*X/nu) * V/U_inf");
xlim([eta(1), eta(end)]);
legend(titles, 'location', 'best');
hold off;

%% Blasius Skin Friction

rho = 1.189;
mu = rho * nu;
REx = (X_E(21:end) - 0.1) * U_inf / nu;

Cf_calc = 0.664 ./ sqrt(REx);
Cf_sim = mu * DUDY(21:end,1) / (0.5 * rho * U_inf^2);

figure(5);
plot(REx, Cf_calc);
hold on;
plot(REx, Cf_sim);
legend("Calculated", "Simulated");
xlabel("Rex");
ylabel("Cf");
hold off;

%% Boundary Layer thickness

delta_sim = Y_C(1) * ones(NX, 1);
for i=1:NX-1
    delta_sim(i) = Y_C(find(U1(i,:) >= (0.99 * U_inf), 1));
end
delta_sim = delta_sim(21:end) ./ (X_E(21:end)- 0.1)';
delta_calc = 5 ./ sqrt(REx);

delta_p_sim = sum(repmat(filter([1, -1], 1, Y_C), NX-20, 1) .*...
    (1-U1(21:end,:)/U_inf), 2) ./ (X_E(21:end)- 0.1)';
delta_p_calc = 1.7208 ./ sqrt(REx);

theta_sim = sum(repmat(filter([1, -1], 1, Y_C), NX-20, 1) .*...
    (U1(21:end,:)/U_inf) .* (1-U1(21:end,:)/U_inf), 2) ./ (X_E(21:end)- 0.1)';
theta_calc = 0.664 ./ sqrt(REx);

figure(6);
plot(REx, delta_calc);
hold on;
plot(REx, delta_sim);
legend("Calculated", "Simulated");
xlabel("Rex");
ylabel("delta/x");
ylim([0 max(delta_calc)]);
hold off;

figure(7);
plot(REx, delta_p_calc);
hold on;
plot(REx, delta_p_sim);
legend("Calculated", "Simulated");
xlabel("Rex");
ylabel("delta*/x");
ylim([0 max(delta_p_calc)]);
hold off;

figure(8);
plot(REx, theta_calc);
hold on;
plot(REx, theta_sim);
legend("Calculated", "Simulated");
xlabel("Rex");
ylabel("theta/x");
ylim([0 max(theta_calc)]);
hold off;