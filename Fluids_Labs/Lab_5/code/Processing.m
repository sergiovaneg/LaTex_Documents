%% Clear environment

clearvars;
close all;
clc;

%% Channel Bulk Velocity + Reynolds

% Q -> 35 L/s
Q = 35E-3;
B = 0.5;
h = 0.42;
Ub = Q / (B*h);

% Reynolds Number
D = 0.06;
U_inf = Ub;
mu = 1e-3;
rho = 9.98E2;
nu = mu / rho;
Re = D * U_inf / nu;

%% Data acquisition + Domain transformation

load("PSVdata.mat");
res = 3040;

x0_px = 465;
y0_px = 290;

x0_m = x0_px / res;
y0_m = y0_px / res;

Grid_x0m = (x0_px - Grid_Xpx) / res;
Grid_y0m = (y0_px - Grid_Ypx) / res;

D = 0.06;
Grid_x0_D = Grid_x0m / D;
Grid_y0_D = Grid_y0m / D;

%% Pixel displacement -> Physical Velocity

Sx0m = - SXpx / res;
Sy0m = - SYpx / res;

fs = 50;
u0 = Sx0m * fs;
v0 = Sy0m * fs;

u0_Ub = u0 / Ub;
v0_Ub = v0 / Ub;

%% Reynolds Average

U0_Ub = mean(u0_Ub, 3, "omitnan");
V0_Ub = mean(v0_Ub, 3, "omitnan");

%% Visualization

figure(1);
quiver(Grid_x0_D, Grid_y0_D, U0_Ub, V0_Ub);
xlabel("x0 [D]");
ylabel("y0 [D]");
hold on;
rectangle('Position', [-1 -0.5 1 1], 'Curvature',[1 1], 'EdgeColor', 'r');
hold off;

figure(2);
streamslice(Grid_x0_D, Grid_y0_D, U0_Ub, V0_Ub);
xlabel("x0 [D]");
ylabel("y0 [D]");
hold on;
rectangle('Position', [-1 -0.5 1 1], 'Curvature',[1 1], 'EdgeColor', 'r');
hold off;

x_custom = linspace(0, 2.5, 1000);
y_custom = linspace(-1, 1, 1000);
U0_c = interp2(Grid_x0_D, Grid_y0_D, U0_Ub, x_custom, 0);
U0_y1 = interp2(Grid_x0_D, Grid_y0_D, U0_Ub, 0.5, y_custom);
U0_y2 = interp2(Grid_x0_D, Grid_y0_D, U0_Ub, 1.5, y_custom);

figure(3);
plot(x_custom, U0_c);
xlabel("x0 [D]");
ylabel("u0(y=0D) [normalized]");

figure(4);
plot(y_custom, U0_y1);
xlabel("y0 [D]");
ylabel("u0 [normalized]");
hold on;
plot(y_custom, U0_y2);
legend(["x = 0.5D"; "x = 1.5D"]);
hold off;

%% Time/Frequency analysis

figure(5);
r_i = [find(Grid_y0_D(:,1) >= 0, 1, "last");
        find(Grid_x0_D(1,:) >= 0, 1, "last")];
t = 0:1/fs:(length(v0)-1)/fs;
v0_Ub_centerline = reshape(v0_Ub(r_i(1), r_i(2), :), size(t));
plot(t, v0_Ub_centerline);
xlabel("t [s]");
ylabel("v0 [normalized]");

figure(6);
[freq, v0_Ub_fft] = fft_of_v0_Ub_velocity(v0_Ub_centerline, fs);
plot(freq, v0_Ub_fft);
f = freq(v0_Ub_fft == max(v0_Ub_fft));
xlabel("f [Hz]");
ylabel("Frequency Spectrum of v0(0,0)");

Sr = f * D / Ub;

%% Spectrum Overlap

figure(7);
for i=3:size(Grid_y0_D,1)-2
    for j=3:r_i(2)
        v0_Ub_single = reshape(v0(i, j, :), size(t)) / Ub;
        [freq, v0_Ub_fft] = fft_of_v0_Ub_velocity(v0_Ub_single, fs);
        plot(freq, v0_Ub_fft);
        hold on;
    end
end
xlabel("f [Hz]");
ylabel("Frequency Spectrum of v0");

%% Animation

window = 0.25/f;
k = round(window * fs);
N = length(t);

u0_Ub_f = movmean(u0 / Ub, k, 3, "omitnan");
v0_Ub_f = movmean(v0 / Ub, k, 3, "omitnan");

figure(8);
quiver(Grid_x0_D, Grid_y0_D, u0_Ub_f(:,:,1), v0_Ub_f(:,:,1));
xlabel("x0 [D]");
ylabel("y0 [D]");
xlim([0 2.5]);
ylim([-1 1]);
hold on;
streamslice(Grid_x0_D, Grid_y0_D, u0_Ub_f(:,:,1), v0_Ub_f(:,:,1));
legend(["Velocity Field";"Streamslices"]);
hold off;

figure(9)
[curlz, ~] = curl(Grid_x0_D, Grid_y0_D,...
                    u0_Ub_f(:,:,1), v0_Ub_f(:,:,1));
surf(Grid_x0_D, Grid_y0_D, curlz);
view(0,90);
xlabel("x0 [D]");
ylabel("y0 [D]");
zlabel("Curl");
xlim([0 2.5]);
ylim([-1 1]);

pause(1/fs);

for n=2:N
    figure(8);
    quiver(Grid_x0_D, Grid_y0_D, u0_Ub_f(:,:,n), v0_Ub_f(:,:,n));
    xlim([0 2.5]);
    ylim([-1 1]);
    hold on;
    streamslice(Grid_x0_D, Grid_y0_D, u0_Ub_f(:,:,n), v0_Ub_f(:,:,n));
    legend(["Velocity Field";"Streamslices"]);
    hold off;
    
    figure(9)
    [curlz, ~] = curl(Grid_x0_D, Grid_y0_D,...
                        u0_Ub_f(:,:,n), v0_Ub_f(:,:,n));
    surf(Grid_x0_D, Grid_y0_D, curlz);
    view(0,90)
    xlabel("x0 [D]");
    ylabel("y0 [D]");
    zlabel("Curl");
    xlim([0 2.5]);
    ylim([-1 1]);
    zlim([-5 5]);
    
    pause(1/fs);
end