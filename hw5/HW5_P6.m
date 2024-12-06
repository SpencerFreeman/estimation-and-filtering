%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/20/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 6 of problem set 5 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW5-P6')


%% P3
kf_example02a % bring in data

nx = length(xhat0);
n = length(thist) + 1;

nmc = 50; % monte carlo's

[P_obs_10, xtil_mu_10, P_obs_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)

nmc = 1000; % monte carlo's

[P_obs_10, xtil_mu_10, P_obs_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)


% temp = nan(2, 2, nmc);
% for i = 1:nmc
%     temp(:, :, i) = xtil(:, 10, i) * xtil(:, 10, i)';
% end
% P_obs2 = mean(temp, 3);
% 
% temp = nan(2, 2, nmc);
% for i = 1:nmc
%     temp(:, :, i) = xtil(:, 35, i) * xtil(:, 35, i)';
% end
% P_obs2_35 = mean(temp, 3);

% 
% te =xtil_10 * xtil_10';


% %% plotting
% close all
% 
% h = figure;
% h.WindowStyle = 'Docked';
% 
% subplot(2, 1, 1)
% plot(ts, xhats(1, :), 'r*'); hold on
% plot(ts, sqrt(phats(1, :)) .* [1; -1], 'bo')
% plot(ts, xtruehist(:, 1), 'g')
% grid on
% legend('Estimate', '1\sigma', '', 'Truth')
% title('Filter Output')
% ylabel('xhat_1')
% 
% subplot(2, 1, 2)
% plot(ts, xhats(2, :), 'r*'); hold on
% plot(ts, sqrt(phats(4, :)) .* [1; -1], 'bo')
% plot(ts, xtruehist(:, 2), 'g')
% grid on
% ylabel('xhat_2')
% xlabel('Time (s)')
% 
% fprintf('P3\n\txhat(50): %f %f\n\tP(50): %f %f %f %f\nP4\n\tPss: %f %f %f %f\n\tStable Error Dynamics: %s\n', ...
%     xhats(:, 50), phats(:, 50), phat_ss, mat2str(is_stable))





function [P_obs_10, xtil_mu_10, P_obs_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)

xtil = nan(nx, n, nmc); % error for MC's

for j = 1:nmc

[xtruehist,zhist] = kf_truthmodel(Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n - 1);

 
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

end % for

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector

xtil(:, :, j) = xtruehist' - xhats;

end % for


xtil_10 = squeeze(xtil(:, 10, :))';
xtil_mu_10 = mean(xtil_10, 1);
P_obs_10 = cov(xtil_10);

P_10 = reshape(phats(:, 10), size(P0));

xtil_35 = squeeze(xtil(:, 35, :))';
xtil_mu_35 = mean(xtil_35, 1);
P_obs_35 = cov(xtil_35);

P_35 = reshape(phats(:, 35), size(P0));

end % function































