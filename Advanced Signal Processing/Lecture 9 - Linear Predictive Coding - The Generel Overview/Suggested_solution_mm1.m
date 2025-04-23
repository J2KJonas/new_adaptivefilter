clear

% Load the data file
load opgave.dat

% Move 128 samples into the vactor "data". It makes a difference where in
% signal the N samples long segment is picked out.
N=128;
data(1:N) = opgave(500:499+N);

% Create a N-point Hamming Window
ham(1:N) = hamming(N);

% Multiply the data segment with the window
data_h(1:N) = data(1:N).*ham(1:N);

% Model order p
p=12;

% Calculate the first p+1 autocorrelation coefficients
% Do it the manual way
for k=1:p+1
    sum = 0;
    for n=1:(N-k+1)
        sum = sum + (data_h(n) * data_h(n+(k-1)));
    end;
    r_1(k) = sum;
end;

% ...or do it the "smart" way
% auto = xcorr(data_h,data_h);
% Here "auto" is the autocorrelation r(k) for k=[-127;127], i.e., 2*128-1
% elements. We need only 11 coefficients, k=1:11
% r_2(1:(p-1)) = auto(128:128+(p-2));


% Next apply the autocorrelation coefficients in the set of linial equ.
% Ra=-r with model order p.

% First create the autocorrelation vector. Here we use the r_1 coeefcients
r(1:p)=r_1(2:p+1);

% Next create the Toeplitz matrix
R = toeplitz(r_1(1:p));

% Solve the equation system
inv_R = inv(R);
a = -inv_R*r';

% Calculate the parametric spectrum
[h,w] = freqz(1,[1 a'], 1000);
spec = abs(h);
plot(w,spec)

% Calculate the residual variance
sum = 0;
for k=1:p;
    sum = sum + a(k) * r_1(k+1);
end;
E_min = r_1(1) + sum



