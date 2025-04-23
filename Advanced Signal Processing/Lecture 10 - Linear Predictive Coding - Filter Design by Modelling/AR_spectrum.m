clear

format long


% Define the model order p
p=32;

% Load the data file
load Amplitude.dat
% ...and move it to tha variable amp
amp = Amplitude;

% Make a frequency scale from 0 to pi
delta_freq = pi/1024;
freq_axis(1)=0;
for i=2:1024
    freq_axis(i) = freq_axis(i-1) + delta_freq;
end

% Make a frequency scale from 0 to 2*pi
delta_freq = (2*pi)/2048;
freq_axis_2pi(1)=0;
for i=2:2048
    freq_axis_2pi(i) = freq_axis_2pi(i-1) + delta_freq;
end

% First we mirror the spectrum around pi
for i=0:1023
    amp(1025+i) = Amplitude(1024-i);
end

%plot(freq_axis_2pi,amp)
%axis([0 7 0 1.4]);


% Next, calculate the power spectrum
power_spec = amp.*amp;
%plot(freq_axis_2pi,power_spec)
%axis([0 7 0 1.6]);

% The autocorrelation coefficients are calculated from the power spectrum 
% as the inverse DTFT. Since the power-spectrum is not an analytical
% function we need to use instead the IDFT.
%
% Be aware that the power spectrum is symmetric and zero-phase, i.e., it is
% "conjugate symmetric". For that reason, calculating the IDFT should
% result in a time series (i.e., autocorr.-coefficients) which is all real.
% However, it seems that Matlab is not able to figure out that the result
% should be real -- it generates a complex result. And to make things worse,
% some of these complex numbers have non-zero imaginary parts -- probably
% due to internal precision used in the ifft calculation.
%
% Aparently, Matlab is aware of this problem... When calculating the ifft,
% it is therefore possible to specify that it should consider the spectrum
% as conjugate symmetric, and thus the output from the ifft becomes all real.
%
% To do this, use the parameter 'symmetric' in the ifft...

auto_corr = ifft(power_spec,'symmetric');
%stem((0:24), auto_corr(1:25))

% Now, generate the autocorrelation matrix/vector
% Remember that the matrix is Toeplitz
R_p = toeplitz(auto_corr(1:p));
r_p = auto_corr(2:p+1);

% Next, solve the equation systems
a_p = -inv(R_p)*r_p;

% Calculate the gain factor which ensures that the resulting spectrum has
% a 0dB DC gain
G_p = 1;
for i=1:p
    G_p = G_p + a_p(i);
end

% Calculate and plot the desired and the approximate amplitude spectrum
[h_p,w] = freqz(1,[1 a_p'], 1024);

plot(freq_axis,amp(1:1024),freq_axis,G_p*abs(h_p))

