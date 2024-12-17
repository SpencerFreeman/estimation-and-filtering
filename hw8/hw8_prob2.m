clear;clc;close all





%%

% high SNR 
load('case01_data.mat')

n = size(zhist, 1) + 1;
nx = length(xhat0);
nz = size(zhist, 2);
thist = 0:delt:(delt*(n - 1));

% UKF
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_ukf(delt,zhist,xhat0,P0,qC,qA,sigma);

% EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))


%% low SNR 
load('case02_data.mat')

n = size(zhist, 1) + 1;
nx = length(xhat0);
thist = 0:delt:(delt*(n - 1));

% xhat a -- UKF
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_ukf(delt,zhist,xhat0_a,P0,qC,qA,sigma);

% xhat a -- EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_a,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

% xhat b -- UKF
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_ukf(delt,zhist,xhat0_b,P0,qC,qA,sigma);

% xhat b -- EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_b,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))


%% plotting
close all

% time histories
names = ["Phase perturbation(rad)", ...
    "Doppler (rad/s)", ...
    "Doppler Rate (rad/s^2)", ...
    "Amplitude"];

f = figure;
f.WindowStyle = 'Docked';
for i = 1:nx
    subplot(nx, 1, i)
    % plot(ts, xhats(i, :), '*'); hold on; grid on
    plot(thist, xhathist(:, i), 'o'); hold on; grid on
    ylabel(names(i))
    if i == 1
        title('Estimated State Time Histories')
        % legend('EKF', 'UKF')
    end % if
end % for
xlabel('Time (s)')

grid on
























