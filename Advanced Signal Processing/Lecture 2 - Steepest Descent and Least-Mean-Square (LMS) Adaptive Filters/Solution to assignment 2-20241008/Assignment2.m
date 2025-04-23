% Assignment 2: NLMS
% Author: Iván López-Espejo (ivl@es.aau.dk)

clear; close all; clc

%% Load Signals
% 'local.asc'   : Clean signal (desired output)
% 'remota.asc'  : Primary input (reference signal with echo)
% 'signal.asc'  : Secondary signal (received signal with echo + desired)

load local.asc       % Clean signal
load remota.asc      % Reference signal (echo source)   Noise
remota = remota';    % Transpose to column vector
load signal.asc      % Primary signal (contains echo)    Signal + noise
%sound(local, 8000)

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

%% Step 3a: NLMS with Double-Talk Detection (DTD)
% Prevent filter from updating after sample 2200 (simulated double-talk)
SNRs = zeros(length(mus), length(ps));
for i = 1:length(mus)
    for j = 1:length(ps)
        mu = mus(i);
        p = ps(j);

        w = zeros(p,1);
        x = zeros(p,1);
        e = zeros(size(remota));

        for k = 1:length(remota)
            x(2:end) = x(1:end-1);
            x(1) = remota(k);
            sig = sum(x.^2);
            e(k) = signal(k) - w' * x;
            if k <= 2200
                w = w + (mu / (sig + 1e-10)) * x * e(k); % DTD active
            end
        end
        SNRs(i,j) = 10 * log10(sum(local.^2) / sum((local - e).^2));
    end
end

%% Step 3b: Plot Weights for Best DTD Configuration
maximum = max(SNRs(:));
[posmu, posp] = find(SNRs == maximum);
opt_mu = mus(posmu);
opt_p = ps(posp);

w = zeros(opt_p, length(remota)+1);
x = zeros(opt_p, 1);
e = zeros(size(remota));

for k = 1:length(remota)
    x(2:end) = x(1:end-1);
    x(1) = remota(k);
    sig = sum(x.^2);
    e(k) = signal(k) - w(:,k)' * x;
    if k <= 2200
        w(:,k+1) = w(:,k) + (opt_mu / (sig + 1e-10)) * x * e(k);
    else
        w(:,k+1) = w(:,k);  % Freeze weights during double-talk
    end
end

% Plot Weights
figure
plot(w')
grid on
title('Weights (NLMS with DTD)')
xlabel('Cycle (n)')
ylabel('Magnitude')
legend(arrayfun(@(i) sprintf('w_%d(n)', i-1), 1:opt_p, 'UniformOutput', false))
disp(['Best SNR (With DTD): ', num2str(maximum), ' dB'])

% % Plot Filtered Output vs Original
% figure
% plot(signal)
% hold on
% plot(e, 'r--')
% hold off
% legend('Received', 'Filtered (DTD)')
% grid on
% xlabel('Time (n)')
% ylabel('Amplitude')
% title('Echo Cancellation with DTD')

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
