% Assignment 1: Wiener filtering.

clear;

nsampl = 10000;  % We generate 10000 samples of the process.
v = randn(nsampl+1,1);
e = randn(nsampl+1,1);
u = zeros(nsampl+1,1);
y = zeros(nsampl+1,1);
for i = 2:nsampl+1
    u(i) = 0.75*u(i-1) + v(i);
    y(i) = u(i) + e(i) + 0.5*e(i-1);
end
u = u(2:end);
y = y(2:end);

% We calculate correlation coefficients.
ruy = xcorr(u,y);
ruy_0 = ruy(nsampl);
ruy_1 = ruy(nsampl+1);
ryy = xcorr(y);
ryy_0 = ryy(nsampl);
ryy_1 = ryy(nsampl+1);

% Weights.
R = toeplitz([ryy_0 ryy_1]);
p = [ruy_0; ruy_1];
w = R \ p;

% We filter the signal.
u_est = filter(w,1,y);

% We plot the result (First 100 samples only).
figure
plot(y(1:100))
hold on
plot(u(1:100),'r')
plot(u_est(1:100),'k--')
hold off
grid on
legend('Noisy input','Desired','Estimated')

% And error variances.
mean((u-y).^2)
mean((u-u_est).^2)