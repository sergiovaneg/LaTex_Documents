%% 0) Env Init

close all;
clearvars;
clc;

%% 1) Channel Bulk Velocity and Reynolds Number

Q = 75E-3;
B = 0.5;
h = 0.45;
Ub = Q/(B*h);

D = 6E-2;
U_inf = Ub;
rho = 998;
mu = 1E-3;
nu = mu/rho;
Re = D*U_inf/nu;

fprintf("Bulk velocity: %d m/s. Reynolds number: %d.\n", Ub, Re);

%% 2) Buoyancy Force

load("FORCEdata_stillwater.mat");
Fx_s = cumsum(Fx)./(1:length(Fx));
Fy_s = cumsum(Fy)./(1:length(Fy));

FB = mean(Fy_s(round(length(Fy_s)/2):end));

fprintf("Buoyancy force: %d N.\n", FB);

%% 3) Reynolds' Average

load("FORCEdata_flow.mat");
FD_avg = cumsum(Fx) ./ (1:length(Fx));
FD_avg = mean(FD_avg(round(length(FD_avg)/4):end));
FL_avg = cumsum(Fy-FB) ./ (1:length(Fy));
FL_avg = mean(FL_avg(round(length(FL_avg)/4):end));

L = 0.185;
A = L * D;
CD = FD_avg / (0.5*rho*A*U_inf^2);
CL = FL_avg / (0.5*rho*A*U_inf^2);

fprintf("Reynolds-averaged Drag force: %d N.\n", FD_avg);
fprintf("Reynolds-averaged Lift force: %d N.\n", FL_avg);
fprintf("Reynolds-averaged Drag coefficient: %d.\n", CD);
fprintf("Reynolds-averaged Lift coefficient: %d.\n", CL);

%% 4) Uncertainties

uU = sqrt(((1/(B*h))*1e-2*Q)^2 + ((-Q/(B^2*h))*5e-4)^2 + ((-Q/(B*h^2))*1e-3)^2);
uA = sqrt((D * 5e-3)^2 + (L * 1e-3)^2);

uD = sqrt(((1/(0.5*rho*A*U_inf^2)) * 45E-3)^2 ...
    + ((-FD_avg/(0.5*rho*A^2*U_inf^2))*uA)^2 ...
    + ((-2*FD_avg/(0.5*rho*A*U_inf^3))*uU)^2);

uL = sqrt(((1/(0.5*rho*A*U_inf^2)) * 25E-3)^2 ...
    + ((-FL_avg/(0.5*rho*A^2*U_inf^2))*uA)^2 ...
    + ((-2*FL_avg/(0.5*rho*A*U_inf^3))*uU)^2);

fprintf("Reynolds-averaged Drag coefficient uncertainty: %d.\n", uD);
fprintf("Reynolds-averaged Lift coefficient uncertainty: %d.\n", uL);

%% 5) Vortex Shedding Frequency

fs = 200;
[freq,FL_fft] = fft_of_force(Fy-FB,fs);
f = freq(FL_fft == max(FL_fft));

fprintf("Vortex-shedding frequency: %d Hz.\n", f);

%% 6) System Natural Frequency

load("FORCEdata_NatOsc.mat");
[freq_x,Fx_fft] = fft_of_force(Dati(1,:),fs);
f_nat_x = freq_x(Fx_fft == max(Fx_fft));
[freq_y,Fy_fft] = fft_of_force(Dati(2,:)-FB,fs);
f_nat_y = freq_y(Fy_fft == max(Fy_fft));
f_nat = mean([f_nat_x, f_nat_y]);

fprintf("Natural frequency of the balance structure: %d Hz.\n", f_nat);

%% 7) Lift-Force Amplitude

FL_filt = filter_force_signal(Fy-FB,fs,2*f);
t=(1:length(Fy))/fs;
figure(1);
plot(t, Fy-FB);
hold on;
plot(t, FL_filt, 'LineWidth', 2);
xlabel("Time [s]");
ylabel("Force [N]");
xlim([t(end)-10,t(end)]);
ylim([-1.5,1.5]);
legend("Original Signal", "Filtered Signal");
hold off;

aL_peak = (mean(findpeaks(FL_filt)) + mean(findpeaks(-FL_filt))) / 2;
aL_rms = sqrt(2) * rms(FL_filt);

fprintf("Lift amplitude from averaged peaks: %d N.\n", aL_peak);
fprintf("Lift amplitude from RMS: %d N.\n", aL_rms);

CL_p_peak = aL_peak/(0.5*rho*A*U_inf^2);
CL_p_rms = aL_rms/(0.5*rho*A*U_inf^2);

fprintf("Fluctuating Lift coefficient from averaged peaks: %d.\n", CL_p_peak);
fprintf("Fluctuating Lift coefficient from RMS: %d.\n", CL_p_rms);

%% 8) Water level error estimate

dW = 2*(2E-3*180E-3)*2E-3;
g = 9.80665;
dFB = dW * rho * g;

dCD = 0;
dCL = -dFB / (0.5*rho*A*U_inf^2);
dCL_p = (sqrt(2) * rms(FL_filt-dFB) - aL_rms)/(0.5*rho*A*U_inf^2);

fprintf("Reynolds-averaged Drag coefficient variation: %d.\n", dCD);
fprintf("Reynolds-averaged Lift coefficient variation: %d.\n", dCL);
fprintf("Fluctuating Lift coefficient variation: %d.\n", dCL_p);

error_CD = abs(dCD/CD);
error_CL = abs(dCL/CL);
error_CL_p = abs(dCL_p/CL_p_rms);

fprintf("Reynolds-averaged Drag coefficient relative error: %d.\n", error_CD);
fprintf("Reynolds-averaged Lift coefficient relative error: %d.\n", error_CL);
fprintf("Fluctuating Lift coefficient relative error: %d.\n", error_CL_p);