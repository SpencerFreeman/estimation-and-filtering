%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/04/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 3 of problem set 5 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW5-P3')


%%
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
    pbar = Fk * phat * Fk' + Qk*1e-0; % propagate state covariance

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


%% plotting
close all

h = figure;
h.WindowStyle = 'Docked';

subplot(2, 1, 1)
plot(ts, xhats(1, :), 'r*'); hold on
plot(ts, sqrt(phats(1, :)) .* [1; -1], 'bo')
grid on
legend('Estimate', '1-\sigma')
title('Filter Output')
ylabel('xhat_1')

subplot(2, 1, 2)
plot(ts, xhats(2, :), 'r*'); hold on
plot(ts, sqrt(phats(4, :)) .* [1; -1], 'bo')
grid on
ylabel('xhat_2')
xlabel('Time (s)')

fprintf('\nxhat(50): %f %f\nP(50): %f %f %f %f\n', xhats(:, 50), phats(:, 50))





































