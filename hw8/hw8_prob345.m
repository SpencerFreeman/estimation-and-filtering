clear;clc;close all


load('measdata_pfexample.mat')

n = length(zkhist); % samples
nx = length(xhat0);
nv = size(Q, 1);
thist = 1:n;

%%
[xhathist_pft,Phist,sigmahist,enuhist] = ...
    particle_filter(zkhist, xhat0, P0, Q, R);

[xhathist_ekf,Phist,sigmahist,enuhist] = ...
    efk(zkhist, xhat0, P0, Q, R);

[xhathist_ukf,Phist,sigmahist,enuhist] = ...
    class_example_ukf(zkhist, xhat0, P0, Q, R);


%% plotting
close all

% time histories
names = ["x1"];

fig = figure;
fig.WindowStyle = 'Docked';
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist_ekf(:, i), '*'); hold on; grid on
    plot(thist, xhathist_pft(:, i), '*'); hold on; grid on
    plot(thist, xhathist_ukf(:, i), '*'); hold on; grid on
    ylabel(names(i))
    if i == 1
        title('Estimated State Time Histories')
        legend('EKF', 'Partical Filter', 'UKF')
    end % if
end % for
xlabel('Time (s)')
grid on













