%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/04/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 3 and 4 of problem set 5 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW5-P3-P4')


%% P3
kf_example02a % bring in data

n = length(thist) + 1;
nx = length(xhat0);
 
t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance

ts = nan(1, n);
xhats = nan(nx, n);
phats = nan(nx * nx, n);
for i = 1:(n - 1)

    ts(i) = t;
    xhats(:, i) = xhat;
    phats(:, i) = phat(:); % unwrap to column vector

    t = thist(i); % s
    xbar = Fk * xhat; % propagate state estimate
    pbar = Fk * phat * Fk' + Gammak * Qk * Gammak'; % propagate state covariance

    zbar = Hk * xbar; % expected measurement 
    z = zhist(i); % actual measurement
    v = z - zbar; % filter innovation

    S = Hk * pbar * Hk' + Rk; % expected measurement covariance
    W = pbar * Hk' * inv(S); % filter gain

    xhat = xbar + W * v; % updated state estimate
    phat = pbar - W * S * W'; % updates state covariance

end

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector


%% P4

[W_ss, pbar_ss, phat_ss] = dlqe(Fk, Gammak, Hk, Qk, Rk);

Fk_error = (eye(size(Fk)) - W_ss*Hk)*Fk;

error_eig = eig(Fk_error);

is_stable = all(abs(error_eig) < 1); % eigenvalues complex magnitudes are stable


%% plotting
close all

h = figure;
h.WindowStyle = 'Docked';

subplot(2, 1, 1)
plot(ts, xhats(1, :), 'r*'); hold on
plot(ts, sqrt(phats(1, :)) .* [1; -1], 'bo')
plot(ts([1, end]), sqrt(phat_ss(1)) .* [1, -1; 1, -1], 'm')
grid on
legend('Estimate', '1\sigma', '', '1\sigma (ss)')
title('Filter Output')
ylabel('xhat_1')

subplot(2, 1, 2)
plot(ts, xhats(2, :), 'r*'); hold on
plot(ts, sqrt(phats(4, :)) .* [1; -1], 'bo')
plot(ts([1, end]), sqrt(phat_ss(4)) .* [1, -1; 1, -1], 'm')
grid on
ylabel('xhat_2')
xlabel('Time (s)')

fprintf('P3\n\txhat(50): %f %f\n\tP(50): %f %f %f %f\nP4\n\tPss: %f %f %f %f\n\tStable Error Dynamics: %s\n', ...
    xhats(:, 50), phats(:, 50), phat_ss, mat2str(is_stable))





































