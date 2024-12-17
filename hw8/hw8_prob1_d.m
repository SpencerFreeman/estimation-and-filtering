clear;clc;close all


nz = 2;




%% d -------------------------------------------------

% high SNR 
load('case01_data.mat')

n = size(zhist, 1) + 1;
nx = length(xhat0);
thist = 0:delt:(delt*(n - 1));

% EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

% IEKF
Niter = 4;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

%% low SNR
load('case02_data.mat')

% xhat0 a -- EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_a,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

% xhat0 a -- IEKF
Niter = 4;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_a,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

% xhat0 b -- EKF
Niter = 1;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_b,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))

% xhat0 b -- IEKF
Niter = 4;
[xhathist,Phist,sigmahist,enuhist] = ...
    pll_iekf(delt,zhist,xhat0_b,P0,qC,qA,sigma,Niter);
ev_bar = mean(enuhist(2:end))














