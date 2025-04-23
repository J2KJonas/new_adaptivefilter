clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal parameters and statistics
M   = 30                     ; % length of the filter
wo  = fir1(M-1,0.2)'         ; % optimal filter coefficients
sigv= 0.01                   ; % observation noise variance
sigu= 1                      ; % variance of the white input signal
Ru  = sigu*eye(M)            ; % correlation matrix
rud = Ru*wo                  ; % cross-correlation vector
sigd= sigv+wo'*Ru*wo         ; % power of the desired signal
Jwo = sigv                   ; % minimum of the cost function
% simulation parameters
N   = 10000                  ; % number of iterations
% parameters of the adaptive filter
w   = zeros(M,N+1)           ; % estimated filter coefficients
e   = zeros(1,N)             ; % errors
y   = zeros(1,N)             ; % output of the filter
Leff= 500                    ; % effective window length
lam = 1-1/Leff               ; % forgetting factor
del = sigu                   ; % variance of initial correlation matrix
P   = eye(M)/del             ; % P(0)
% generate the signals
u   = sqrt(sigu)*randn(1,N+M-1); % white gaussian noise
U   = flipud(hankel(u(1:M),u(M:N+M-1))); %matrix of input vectors
z   = filter(wo,1,u)         ; % create noiseless desired signal
d   = z(M:N+M-1)+sqrt(sigv)*randn(1,N); % desired vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the adaptive Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:N
  z    = P*U(:,n)            ; % temporary variable
  k    = z/(lam+U(:,n)'*z)   ; % gain vector
  y(n) = U(:,n)'*w(:,n)      ; % output of the filter
  e(n) = d(n)-y(n)           ; % a priori error
  w(:,n+1) = w(:,n)+k*e(n)   ; % update the filter coefficients
  P    = (P-k*z')/lam        ; % update inverse correlation matrix
end
% these plots are inspired by the adapfilt.rls entry in the matlab manual
subplot(3,1,1);
plot(1:N,[d;y;e]);
legend('Desired','Output','Error');
title('System Identification of FIR Filter');
xlabel('Time Index'); ylabel('Signal Value');
subplot(2,1,2);
stem([wo, w(:,N+1)]);
legend('Actual','Estimated'); grid on;
xlabel('Coefficient #'); ylabel('Coefficient Value'); 
 