clear;clc;close all


load('measdata_pfexample.mat')

n = length(zkhist); % samples
nx = length(xhat0);
nv = size(Q, 1);
thist = 1:n;

%%
Ns = 10;
[xhathist_10,Phist,sigmahist,enuhist] = ...
    particle_smoother(zkhist, xhat0, P0, Q, R, Ns);

Ns = 100;
[xhathist_100,Phist,sigmahist,enuhist] = ...
    particle_smoother(zkhist, xhat0, P0, Q, R, Ns);

Ns = 1000;
[xhathist_1000,Phist,sigmahist,enuhist] = ...
    particle_smoother(zkhist, xhat0, P0, Q, R, Ns);



%% plotting
close all

% time histories
names = ["x1"];

fig = figure;
fig.WindowStyle = 'Docked';
for i = 1:nx
    subplot(nx, 1, i)
    plot(thist, xhathist_10(:, i), '*'); hold on; grid on
    plot(thist, xhathist_100(:, i), '*'); hold on; grid on
    plot(thist, xhathist_1000(:, i), '*'); hold on; grid on
    % plot(thist, xhathist_ukf(:, i), '*'); hold on; grid on
    ylabel(names(i))
    if i == 1
        title('Estimated State Time Histories')
        legend('10', '100', '1000')
    end % if
end % for
xlabel('Time (s)')
grid on













