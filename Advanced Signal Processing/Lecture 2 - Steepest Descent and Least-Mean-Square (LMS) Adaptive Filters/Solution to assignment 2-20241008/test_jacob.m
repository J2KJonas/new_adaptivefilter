% Assignment 2: NLMS
% Author: Iván López-Espejo (ivl@es.aau.dk)

clear; close all; clc

%% Load Signals
[remota, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\primary.wav");   %noise + clean signal
[signal, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\secondary.wav");
remota = remota';    % Transpose to column vector
[local, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");

%% Step 1: Baseline SNR (no echo cancellation)
SNR = 10 * log10(sum(local.^2) / sum((local - signal).^2));
disp(['Baseline SNR: ', num2str(SNR), ' dB'])

%% Step 2a: Grid Search for Best NLMS Parameters (mu, p)
mus = 1e-4 * (2.^(0:9));   % Step sizes to try
ps = 1:100;                 % Filter lengths to try
SNRs = zeros(length(mus), length(ps));  % To store resulting SNRs

for i = 1:length(mus)
    for j = 1:length(ps)
        % Parameters
        mu = mus(i);
        p = ps(j);

        % Initialization
        w = zeros(p, 1);       % Filter weights
        x = zeros(p, 1);       % Filter input buffer
        e = zeros(size(remota));  % Filter output (error signal)

        % Adaptive Filtering
        for k = 1:length(remota)
            x(2:end) = x(1:end-1);   % Shift buffer
            x(1) = remota(k);        % Insert new sample
            sig = sum(x.^2);         % Norm of input
            e(k) = signal(k) - w' * x;
            w = w + (mu / (sig + 1e-10)) * x * e(k); % NLMS update
        end

        % Compute output SNR
        SNRs(i, j) = 10 * log10(sum(local.^2) / sum((local - e).^2));
    end
end

%% Step 2b: Plot Filter Weights for Best NLMS Configuration
maximum = max(SNRs(:));
[posmu, posp] = find(SNRs == maximum);
opt_mu = mus(posmu);     % Optimal mu
opt_p = ps(posp);        % Optimal filter order

% Re-run filtering using best configuration
w = zeros(opt_p, length(remota)+1);  % Matrix to store weights over time
x = zeros(opt_p, 1);
e = zeros(size(remota));

for k = 1:length(remota)
    x(2:end) = x(1:end-1);
    x(1) = remota(k);
    sig = sum(x.^2);
    e(k) = signal(k) - w(:,k)' * x;
    w(:,k+1) = w(:,k) + (opt_mu / (sig + 1e-10)) * x * e(k);
end

% Plot Weights Over Time
figure
plot(w')
grid on
title('Weights (Optimal NLMS)')
xlabel('Cycle (n)')
ylabel('Magnitude')
legend(arrayfun(@(i) sprintf('w_%d(n)', i-1), 1:opt_p, 'UniformOutput', false))
disp(['Best SNR (No DTD): ', num2str(maximum), ' dB'])

% Plot Filtered Output vs Original
figure
plot(signal)
hold on
plot(e, 'r--')
hold off
legend('Received', 'Filtered')
grid on
xlabel('Time (n)')
ylabel('Amplitude')
title('Echo Cancellation Output')

%% Step 4: Frequency Response of Final Filter
wf = w(:, end);  % Final filter weights
fr = 20 * log10(abs(freqz(wf, 1, 257)));  % Frequency response
x_range = 0:4000/256:4000;

figure
plot(x_range, fr)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency Response of Adaptive Filter (DTD)')
grid on
