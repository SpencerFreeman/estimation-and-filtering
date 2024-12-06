%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/04/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 5 of problem set 5 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW6-P2')


%%

kf_example03a % bring in data

[ts, xhats, phats, evs] = ...
    filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, thist, zhist);

[ts, xhats, phats, evs] = ...
    SRIF(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, thist, zhist);


kf_example03b % bring in data

[ts, xhats, phats, evs] = ...
    filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, thist, zhist);

[ts, xhats, phats, evs] = ...
    SRIF(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, thist, zhist);



%% plotting
close all








function [ts, xhats, phats, evs] = ...
    filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, thist, zhist)

n = length(thist) + 1;
nx = length(xhat0);
 
t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
xhats = nan(nx, n);
phats = nan(nx * nx, n);
evs = nan(1, n);
for i = 1:(n - 1)

    ts(i) = t;
    xhats(:, i) = xhat;
    phats(:, i) = phat(:); % unwrap to column vector
    evs(i) = ev;

    t = thist(i); % s
    xbar = Fk * xhat; % propagate state estimate
    pbar = Fk * phat * Fk' + Gammak * Qk * Gammak'; % propagate state covariance

    zbar = Hk * xbar; % expected measurement 
    z = zhist(i); % actual measurement
    v = z - zbar; % filter innovation

    S = Hk * pbar * Hk' + Rk; % expected measurement covariance
    Sinv = inv(S);
    W = pbar * Hk' * Sinv; % filter gain

    ev = v' * Sinv * v; % estimation error statistic

    xhat = xbar + W * v; % updated state estimate
    phat = pbar - W * S * W'; % updates state covariance

end

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector
evs(n) = ev;

end % function








































