%% Exercise 1: Confirm that Phi is an orthonormal set
% Define T
T = 1;

% Define the basis functions phi1 and phi2
phi1 = @(t) sqrt(2) .* (t >= 0 & t < T/2);
phi2 = @(t) sqrt(2) .* (t >= T/2 & t <= T);

% Compute the inner products to check orthogonality and normalization
phi1_phi1 = integral(@(t) phi1(t) .* phi1(t), 0, T);  % <phi1, phi1>
phi2_phi2 = integral(@(t) phi2(t) .* phi2(t), 0, T);  % <phi2, phi2>
phi1_phi2 = integral(@(t) phi1(t) .* phi2(t), 0, T);  % <phi1, phi2>

fprintf('Exercise 1: Orthogonality and Normalization\n');
fprintf('Inner product <phi1, phi1>: %f\n', phi1_phi1);
fprintf('Inner product <phi2, phi2>: %f\n', phi2_phi2);
fprintf('Inner product <phi1, phi2>: %f\n', phi1_phi2);
fprintf('-----------------------------\n\n');

%% Exercise 2: Check that the two signals x(t) and y(t) are in V
% Define the signals x(t) and y(t)
x = @(t) 10 .* (t >= 0 & t < T/2);
y = @(t) 5 .* (t >= 0 & t < T/2) - 10 .* (t >= T/2 & t <= T);

% Compute the coefficients for x(t) and y(t) in terms of phi1 and phi2
x1 = integral(@(t) x(t) .* phi1(t), 0, T);  % coefficient for phi1 in x(t)
x2 = integral(@(t) x(t) .* phi2(t), 0, T);  % coefficient for phi2 in x(t)
y1 = integral(@(t) y(t) .* phi1(t), 0, T);  % coefficient for phi1 in y(t)
y2 = integral(@(t) y(t) .* phi2(t), 0, T);  % coefficient for phi2 in y(t)

fprintf('Exercise 2: Coefficients of x(t) and y(t) in terms of phi1 and phi2\n');
fprintf('x1: %f, x2: %f\n', x1, x2);
fprintf('y1: %f, y2: %f\n', y1, y2);
fprintf('-----------------------------\n\n');

%% Exercise 3: Example of a signal not in V
% Define a signal z(t) that is not in V
z = @(t) sin(pi * t / T);

% Plot the signal z(t)
t_values = linspace(0, T, 1000);
figure;
plot(t_values, z(t_values));
title('Exercise 3: Signal z(t) = sin(\pi t / T), which is not in V');
xlabel('t');
ylabel('z(t)');
grid on;

fprintf('Exercise 3: Plotted the signal z(t) = sin(pi t / T), which is not in V.\n');
fprintf('-----------------------------\n\n');

%% Exercise 4: Represent x(t) and y(t) in the basis Phi
% The coefficients for x(t) and y(t) have already been computed in Exercise 2.
fprintf('Exercise 4: Representing x(t) and y(t) in the basis Phi\n');
fprintf('x(t) = %f * phi1(t) + %f * phi2(t)\n', x1, x2);
fprintf('y(t) = %f * phi1(t) + %f * phi2(t)\n', y1, y2);
fprintf('-----------------------------\n\n');

%% Exercise 5: Evaluate <x, y> and compare with sum(x_n * y_n)
% Compute the inner product <x, y>
inner_product_xy = integral(@(t) x(t) .* y(t), 0, T);  % <x, y>

% Compute the sum of the products of the coefficients
sum_coefficients = x1 * y1 + x2 * y2;

fprintf('Exercise 5: Inner product <x, y> and sum of coefficient products\n');
fprintf('Inner product <x, y>: %f\n', inner_product_xy);
fprintf('Sum of products of coefficients: %f\n', sum_coefficients);
fprintf('-----------------------------\n\n');

%% Exercise 6: Construct a signal z(t) in V which is orthogonal to x(t)
% To be orthogonal, <z, x> = 0
% Solve for z1 and z2 such that z1 * x1 + z2 * x2 = 0
z1 = 0;  % Set z1 = 0 since it simplifies the problem
z2 = 1;  % Choose any non-zero value for z2

% Define the orthogonal signal z(t)
z_orthogonal = @(t) z1 * phi1(t) + z2 * phi2(t);

% Verify orthogonality <z, x>
inner_product_zx = integral(@(t) z_orthogonal(t) .* x(t), 0, T);  % <z, x>

fprintf('Exercise 6: Constructing a signal z(t) orthogonal to x(t)\n');
fprintf('Inner product <z, x> (should be 0): %f\n', inner_product_zx);
fprintf('-----------------------------\n\n');
