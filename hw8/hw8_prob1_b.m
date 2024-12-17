clear;clc;close all


nz = 1;


%% b -------------------------------------------------

% high SNR 
load('case01_data.mat')

xhat0 = xhat0(1:3);
P0 = P0(1:3, 1:3);

n = size(zhist, 1) + 1;
nx = length(xhat0);
thist = 0:delt:(delt*(n - 1));

[xhathist, Phist, sigmahist, nhist, enuhist] = ...
    pll_lkf(delt, zhist, xhat0, P0, qC, sigma);

ev_bar = mean(enuhist(2:end))

% bup = chi2inv(0.95, nz)
% blo = chi2inv(0.05, nz)

alpha = .05;
bup = chi2inv(1 - (1 - alpha)/2, 2);
blo = chi2inv((1 - alpha)/2, 2);


% figure
% histogram(enuhist)
% grid on

figure
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist(:, i), '.')
    grid on
end % for
xlabel('Time (s)')

% low SNR
load('case02_data.mat')

xhat0_a = xhat0_a(1:3);
xhat0_b = xhat0_b(1:3);
P0 = P0(1:3, 1:3);

% xhat0 a
[xhathist, Phist, sigmahist, nhist, enuhist] = ...
    pll_lkf(delt, zhist, xhat0_a, P0, qC, sigma);

ev_bar = mean(enuhist(2:end))

figure
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist(:, i), '.')
    grid on
end % for
xlabel('Time (s)')


% xhat0 b
[xhathist, Phist, sigmahist, nhist, enuhist] = ...
    pll_lkf(delt, zhist, xhat0_b, P0, qC, sigma);

ev_bar = mean(enuhist(2:end))

figure
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist(:, i), '.')
    grid on
end % for
xlabel('Time (s)')

















