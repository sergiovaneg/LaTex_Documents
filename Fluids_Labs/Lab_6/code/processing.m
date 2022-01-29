%% Env re-init

close all;
clearvars;
clc;

%% Grid-Independence

RANS_files = dir("RANS_*.mat");
Force_files = dir("forces-s_*.csv");

labels = ["3° grid-spacing";
            "2° grid-spacing";
            "1° grid-spacing";
            "Separation point"];
rho = 998.23;
mu = 1.006E-3;
Dc = 0.06;
H = 1;
U_inf = 0.4;
separation = zeros(3,length(RANS_files));

for i=1:length(RANS_files)
    % Import
    load(RANS_files(i).name);
    F = readmatrix(Force_files(i).name);

    figure(1);
    plot(X_C(1:NX/2) * 180 / pi, P1(1:NX/2,2));
    hold on;

    aux = filter([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560],...
                    1, P1(1:NX/2,2));
    aux(1:NX/5) = 1;
    idx = find(aux < 0, 1) - 4;
    separation(1,i) = X_C(idx) * 180 / pi;
    separation(2,i) = P1(idx,2);
    
    aux = filter([1 -1], 1, STRS(1:NX/3,2));
    aux(1:NX/6) = -1;
    idx = find(aux>=-1e-6, 1);
    aux = rho*STRS(1:NX/2,2);
    aux(idx:end) = -aux(idx:end);

    figure(2);
    plot(X_C(1:NX/2) * 180 / pi, rho*STRS(1:NX/2,2));
    hold on;

    m = (aux(idx) - aux(idx-1)) / (X_C(idx) - X_C(idx-1));
    separation(3,i) = (X_C(idx-1) - aux(idx-1) / m) * 180 / pi;

    figure(3);
    plot(F(:,1), (F(:,4) / H) / (0.5 * rho * Dc * U_inf^2));
    hold on;
end

figure(1);
xlim([0,180]);
xlabel("x [°]");
ylabel("Pressure [Pa]");
scatter(separation(1,:), separation(2,:));
legend(labels);
hold off;

figure(2);
xlim([0,180]);
xlabel("x [°]");
ylabel("Wall Shear-stress [Pa]");
scatter(separation(3,:), zeros(1,length(RANS_files)));
legend(labels);
hold off;

figure(3);
xlim([0,2E3]);
xlabel("Iterations [normalized]");
ylabel("Drag coefficient [adimensional]");
legend(labels(1:3));
hold off;

%% Dynamic sensitivity analysis

F_06 = readmatrix("forces_dynamic_06ms.csv");
F_10 = readmatrix("forces_dynamic_10ms.csv");
F_20 = readmatrix("forces_dynamic_20ms.csv");

figure(4);
plot(F_06(:,1),F_06(:,4));
hold on;
plot(F_10(:,1), F_10(:,4));
plot(F_20(:,1), F_20(:,4));
legend("dt=6 ms","dt=10 ms","dt=20 ms");
xlabel("time [s]");
ylabel("Lift Force [N]");
hold off;

figure(5);
semilogy(F_06(:,1), abs(F_06(:,5)));
hold on;
semilogy(F_10(:,1), abs(F_10(:,5)));
semilogy(F_20(:,1), abs(F_20(:,5)));
legend("dt=6 ms","dt=10 ms","dt=20 ms");
xlabel("time [s]");
ylabel("Drag Force [N]");
hold off;

%% Drag/Lift Coefficient Evolution

Fx = F_06(:,4);
Fy = -F_06(:,5);

H = 1;
Cd = (Fy/H) / (0.5*rho*Dc*U_inf^2);
Cl = (Fx/H) / (0.5*rho*Dc*U_inf^2);

figure(6);
subplot(2,1,1);
plot(F_06(:,1), Cd);
xlabel("time [s]");
ylabel("Drag Coefficient []");
subplot(2,1,2);
plot(F_06(:,1), Cl);
xlabel("time [s]");
ylabel("Lift Coefficient []");

%% Time-Averaged results

Fx_avg = mean(Fx(end/4:end));
Fy_avg = mean(Fy(end/4:end));

Cd_avg = (Fy_avg/H) / (0.5*rho*Dc*U_inf^2);
Cl_avg = (Fx_avg/H) / (0.5*rho*Dc*U_inf^2);

%% Strouhal Number

dt = F_06(2,1) - F_06(1,1);
FX = fft(Fx(end/4:end)-Fx_avg);
freq = 0:1/(dt*length(FX)):1/dt-1/(dt*length(FX));
f_wake = freq(find(abs(FX)==max(abs(FX)),1));
St = f_wake*Dc/U_inf;