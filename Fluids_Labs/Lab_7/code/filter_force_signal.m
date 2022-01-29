function [signal_force_filtered] = filter_force_signal(signal_force,fs,fco)

%INPUT
% - signal_force            time history of force values at all time steps
% - fs                      frequency of acquisition of the force signal
% - fco                     cut-off frequency

%OUTPUT
% - signal_force_filtered   time history of force values after applying the filter


ny=fs/2; %%Nyquist frequency
Wn=fco/ny; %%dimensionless cut-off frequency
N=4; %%order of the filter (the higher the order, %%the sharper the frequency cut)
%Butterworth low-pass filter
[b,a]=butter(N,Wn,'low');
signal_force_filtered = filtfilt(b,a,signal_force);    