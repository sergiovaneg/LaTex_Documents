function [freq,force_fft] = fft_of_force(signal_force,fs)

%INPUT
% - signal_force        time history of force values at all time steps
% - fs                  frequency of acquisition of the force signal

%OUTPUT
% - freq                vector of frequencies in the FFT calculation
% - force_fft           FFT of force signal (zero-mean)


dt = 1/fs; %time step (fs: sampling frequency)
ns=numel(signal_force-mean(signal_force));   %Number of samples
time_rec=[0:ns-1]/fs;

% FFT of the velocity signal
freq_fft = linspace(-fs/2,fs/2,length(time_rec));
signal_fft = fftshift(fft(signal_force,length(freq_fft)))./fs;
freq = freq_fft(freq_fft>=0);
force_fft = abs(signal_fft(freq_fft>=0));