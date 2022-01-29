function [freq,v0_Ub_fft] = fft_of_v0_Ub_velocity(signal_v0_Ub_rec,fs)

%INPUT
% - signal_v_rec        time history of v0 values at all time steps (might include NaN entries)
% - fs                  frequency of acquisition of the v0 signal

%OUTPUT
% - freq                vector of frequencies in the FFT calculation
% - v0_fft              FFT of v0 signal


dt = 1/fs; %time step (fs: sampling frequency)
time_rec=dt:dt:dt*length(signal_v0_Ub_rec); 

%Eliminate NaN entries from the signal
indexes_wo_nan = ~isnan(signal_v0_Ub_rec);
time_wo_nan = time_rec(indexes_wo_nan);
signal_v0_Ub_wo_nan = signal_v0_Ub_rec(indexes_wo_nan);

%Reinterpolate signal without NaN entries 
signal_v0_Ub = interp1(time_wo_nan,signal_v0_Ub_wo_nan,time_rec);
signal_v0_Ub = signal_v0_Ub - mean(signal_v0_Ub);

% FFT of the velocity signal
freq_fft = linspace(-fs/2,fs/2,length(time_rec));
signal_fft = fftshift(fft(signal_v0_Ub,length(freq_fft)))./fs;
freq = freq_fft(freq_fft>=0);
v0_Ub_fft = abs(signal_fft(freq_fft>=0));