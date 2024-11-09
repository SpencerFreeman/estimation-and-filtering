%% Implement a Kalman filter for the example problem that was presented in class
% Spencer Freeman, 11/04/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 5 of problem set 5 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW5-P5')


%% P3
kf_example02b % bring in data

Qks(1)     =   40.00000000000000;                         % for all k
Qks(2)     =    0.40000000000000;     % alternate          % for all k
Qks(3)     =    0.00400000000000;

evbars = nan(1, length(Qks));
xhats_save = cell(1, 3);
phats_save = cell(1, 3);
for i = 1:length(Qks)
    [ts, xhats, phats, evs] = ...
        filter(Fk, Gammak, Hk, Qks(i), Rk, xhat0, P0, thist, zhist);

    n = length(evs);
    evbars(i) = mean(evs)*n;
    xhats_save{i} = xhats;
    phats_save{i} = phats;
end

prob = 0.99;

bup = chi2inv(1 - (1 - prob)/2, n);
blo = chi2inv((1 - prob)/2, n);

chi2cdf(bup, n) - chi2cdf(blo, n) % should equal prob (0.99)

i_best = 2;

rms_x1 = rms(xhats_save{i_best}(1, (end - 40):end) - xhats_save{1}(1, (end - 40):end));
rms_v1 = rms(xhats_save{i_best}(2, (end - 40):end) - xhats_save{1}(2, (end - 40):end));

rms_x2 = rms(xhats_save{i_best}(1, (end - 40):end) - xhats_save{3}(1, (end - 40):end));
rms_v2 = rms(xhats_save{i_best}(2, (end - 40):end) - xhats_save{3}(2, (end - 40):end));

sig_x = sqrt(phats_save{i_best}(1, end));
sig_v = sqrt(phats_save{i_best}(4, end));

%% plotting
close all

h1 = figure;
h1.WindowStyle = 'Docked';

x1 = linspace(0, 100, 100);
y1 = chi2pdf(x1, n);
plot(x1, y1); hold on
xline(blo)
xline(bup)
plot(evbars, zeros(size(evbars)), '*')
grid on

% h2 = figure;
% h2.WindowStyle = 'Docked';
% 
% subplot(2, 1, 1)
% plot(ts, xhats(1, :), 'r*'); hold on
% plot(ts, sqrt(phats(1, :)) .* [1; -1], 'bo')
% grid on
% legend('Estimate', '1\sigma', '')
% title('Filter Output')
% ylabel('xhat_1')
% 
% subplot(2, 1, 2)
% plot(ts, xhats(2, :), 'r*'); hold on
% plot(ts, sqrt(phats(4, :)) .* [1; -1], 'bo')
% grid on
% ylabel('xhat_2')
% xlabel('Time (s)')

fprintf('P5\n\tLower Threshold: %f \n\tUpper Threshold: %f \n', blo, bup)
fprintf('\n\tBest Qk: %f \n', Qks(i_best))
fprintf('\n\tRMS x1/sigma_x: %f \n\tRMS x2/sigma_x: %f \n\tRMS v1/sigma_v: %f \n\tRMS v2/sigma_v: %f \n', ...
    rms_x1/sig_x, rms_x2/sig_x, rms_v1/sig_v, rms_v2/sig_v)






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








































