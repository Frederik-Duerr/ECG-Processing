function [filtered_signal] = bandpass_filter(Fs, signal, num_leads)
% BANDPASS_FILTER Apply bandpass and 50 Hz notch filter to multi-lead ECG
% Inputs:
%   Fs        - Sampling frequency [Hz]
%   signal    - ECG signal matrix [samples x leads]
%   num_leads - Number of leads
% Output:
%   filtered_signal - Filtered ECG

%% -----------------------------
% 1. Bandpass Filter (0.5-40 Hz)
%% -----------------------------

f_low  = 0.5;   % Lower cutoff [Hz]
f_high = 40;    % Upper cutoff [Hz]
[b_bp, a_bp] = butter(3, [f_low f_high]/(Fs/2));  % 3rd order Butterworth

%% -----------------------------
% 2. Notch Filter at 50 Hz
%% -----------------------------

wo = 50/(Fs/2);       % Normalized frequency
bw = wo/35;           % Bandwidth (adjustable)
[b_notch, a_notch] = iirnotch(wo, bw);

%% -----------------------------
% 3. Apply filters to each lead
%% -----------------------------

filtered_signal = zeros(size(signal));
for i = 1:num_leads
    % 1. Bandpass
    temp = filtfilt(b_bp, a_bp, signal(:,i));
    % 2. Notch at 50 Hz
    filtered_signal(:,i) = filtfilt(b_notch, a_notch, temp);
end

end
