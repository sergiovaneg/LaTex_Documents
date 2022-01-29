%% Data retrieval

clearvars;
close all;
Fs = 100;
files = dir("calibr*.mat");
forces = [0., -0.314, -0.411, -0.695, -0.891, -1.186, -1.676, -2.656,...
            -5.598, -10.501, 0.126, 1.097, 2.078, 4.039, 8.942, 13.845];
voltages = zeros(1, length(files));
min_deltas = voltages;
figure();
for i = 1:length(files)
    load(files(i).name);
    deltas = filter([1 -1], 1, sort(V));
    deltas = deltas(2:end);
    min_deltas(i) = min(deltas(deltas > 0));
    averaged_data = cumsum(V) ./ (1:length(V));
    loglog((1:length(averaged_data))/Fs, ((averaged_data-averaged_data(end)) / averaged_data(end)).^2);
    hold on;
    voltages(i) =...
        median(averaged_data(10*Fs:end));
end
xlabel('time(s)');
ylabel('Averaged sensor reading (Normalized)');
ylim([1e-6, 1]);
hold off;
figure();
scatter(forces, voltages);

%% Curve fitting

p = polyfit(forces, voltages, 1);
m = p(1);
b = p(2);
voltages_fit = m * forces + b;

hold on;
plot(forces, voltages_fit);
xlabel('Force (N)');
ylabel('Voltage (V)');
legend('Data points', 'Linear fitting');
hold off;

%% Uncertainty - L_{inf} approach

U1 = max(abs((voltages-b)/m - forces));

%% Uncertainty - Confidence Intervals

Sf_sq = sum(((voltages - b) / m - forces).^2) / (length(forces)-2);
U2 = 2 * sqrt(Sf_sq);

%% Resolution

voltage_res = min(min_deltas);
force_res = abs(voltage_res / m);