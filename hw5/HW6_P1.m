%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/20/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 1 of problem set 6 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW6-P1')


%% P3
kf_example02a % bring in data

nx = length(xhat0);
n = length(thist) + 1;

nmc = 50; % monte carlo's

[P_obs_10, P_10, xtil_mu_10, P_obs_35, P_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)

nmc = 1000; % monte carlo's

[P_obs_10, P_10, xtil_mu_10, P_obs_35, P_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)

disp('Errors are higher, estimated covariances are much smaller than is actually the case.')







function [P_obs_10, P_10, xtil_mu_10, P_obs_35, P_35, xtil_mu_35] = run_mc( ...
    thist, Fk,Gammak,Hk,Qk,Rk,xhat0,P0,n, nx, nmc)

xtil = nan(nx, n, nmc); % error for MC's

for j = 1:nmc

[xtruehist,zhist] = kf_truthmodel(Fk,Gammak,Hk,40,Rk,xhat0,P0,n - 1);

 
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
    pbar = Fk * phat * Fk' + Gammak * 0.004 * Gammak'; % propagate state covariance

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































