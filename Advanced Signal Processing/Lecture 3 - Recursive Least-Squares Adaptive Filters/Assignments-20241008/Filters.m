clear all;
clc;

% Read audio files
[u, fs1] = audioread('Noisy_Music.wav'); % Noisy signal    Primary
[d, fs1] = audioread('Noise.wav'); % Desired signal        Secondary
[x, fs1] = audioread('Music.wav'); % Original signal       Clean

N = length(u); % Length of the noisy signal
M = 12; % Filter Order

% *************************************************************************
% LMS Filter
tic
mu = 0.01; % Step size
w = randn(M, 1); % Initialize filter weights
padded_u = [zeros(M-1, 1); u]; % Pad input signal
y = zeros(N, 1); % Initialize output

for n = 1:N
    u_vect = padded_u(n:n+M-1); % Current input vector
    e = d(n) - w' * u_vect; % Error signal
    w = w + mu * e * u_vect; % Update weights
    y(n) = w' * u_vect; % Filtered output
end

filtered_signal_LMS = u - y; % Filtered signal
toc

% *************************************************************************
% NLMS Filter
tic
mu = 1; % Step size
w = randn(M, 1); % Initialize filter weights
padded_u = [zeros(M-1, 1); u]; % Pad input signal
y = zeros(N, 1); % Initialize output
Eps = 0.0001; % Small constant for stability

for n = 1:N
    u_vect = padded_u(n:n+M-1); % Current input vector
    mu1 = mu / (Eps + norm(u_vect)^2); % Normalized step size
    e = d(n) - w' * u_vect; % Error signal
    w = w + mu1 * e * u_vect; % Update weights
    y(n) = w' * u_vect; % Filtered output
end

filtered_signal_NLMS = u - y; % Filtered signal
toc

% *************************************************************************
% RLS Filter
tic
lambda = 1 - 1 / (0.1 * M); % Forgetting factor
delta = 0.01; % Small constant for initialization
P = 1 / delta * eye(M); % Initialize inverse correlation matrix
w = randn(M, 1); % Initialize filter weights
padded_u = [sqrt(delta) * randn(M-1, 1); u]; % Pad input signal
y = zeros(N, 1); % Initialize output

for n = 1:N
    u_vect = padded_u(n:n+M-1); % Current input vector
    PI = P * u_vect; % Intermediate calculation
    gain_k = PI / (lambda + u_vect' * PI); % Gain
    prior_error = d(n) - w' * u_vect; % Error signal
    w = w + prior_error * gain_k; % Update weights
    P = P / lambda - gain_k * (u_vect' * P) / lambda; % Update inverse correlation matrix
    y(n) = w' * u_vect; % Filtered output
end

filtered_signal_RLS = u - y; % Filtered signal
toc

% *************************************************************************
% Visualization of Signals
figure;
subplot(4, 1, 1);
plot(u);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(4, 1, 2);
plot(filtered_signal_LMS);
title('Filtered Signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(4, 1, 3);
plot(filtered_signal_NLMS);
title('Filtered Signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(4, 1, 4);
plot(filtered_signal_RLS);
title('Filtered Signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');

% Save filtered signals to audio files
audiowrite('Filtered_LMS.wav', filtered_signal_LMS, fs1);
audiowrite('Filtered_NLMS.wav', filtered_signal_NLMS, fs1);
audiowrite('Filtered_RLS.wav', filtered_signal_RLS, fs1);
