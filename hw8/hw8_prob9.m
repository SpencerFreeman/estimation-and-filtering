clear;clc;close all





kf_example02a % bring in data

nx = length(xhat0);

%%

[xhathist_lkf,Phist,sigmahist,enuhist] = ...
    prob9_lkf(zhist, xhat0, P0, Qk, Rk, Gammak, Fk, Hk);

Ns = 1;
[xhathist_pft_1,Phist,sigmahist,enuhist] = ...
    prob9_particle_filter(zhist, xhat0, P0, Qk, Rk, Gammak, Fk, Hk, Ns);

Ns = 10;
[xhathist_pft_10,Phist,sigmahist,enuhist] = ...
    prob9_particle_filter(zhist, xhat0, P0, Qk, Rk, Gammak, Fk, Hk, Ns);

Ns = 100;
[xhathist_pft_100,Phist,sigmahist,enuhist] = ...
    prob9_particle_filter(zhist, xhat0, P0, Qk, Rk, Gammak, Fk, Hk, Ns);



%% plotting
close all

% time histories
names = ["x1", "x2"];

fig = figure;
fig.WindowStyle = 'Docked';
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist_lkf(:, i), 'o'); hold on; grid on
    plot(thist, xhathist_pft_1(:, i), '*'); hold on; grid on
    plot(thist, xhathist_pft_10(:, i), '*'); hold on; grid on
    plot(thist, xhathist_pft_100(:, i), '*'); hold on; grid on

    ylabel(names(i))
    if i == 1
        title('Estimated State Time Histories')
        legend('LKF', 'PFT-1', 'PFT-10', 'PFT-100')
    end % if
end % for
xlabel('Time (s)')
grid on























