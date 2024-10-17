clear;clc;close all



load('radarmeasdata_cart.mat')

rhobhist = rhobhist';
rhoahist = rhoahist';
thist = thist';

% x = [ψ0; y1r0; y2r0;ψ;vr].

la = -1; % m
lb = 1; % m

sig_rhoa = .005; % m
sig_rhob = .005; % m

R = diag([sig_rhoa^2, sig_rhob^2]); % m^2
Ra = chol(R); % Cholesky-factor square root, m^2
Rainv = inv(Ra);

y1z = (rhobhist.^2 - rhoahist.^2 - lb^2 + la^2) / 2 / (la - lb); % m
y2z = sqrt(rhoahist.^2 - (y1z - la).^2); % m

v1z = diff(y1z) / (thist(2) - thist(1)); v1z = [v1z, v1z(end)]; % m/s
v2z = diff(y2z) / (thist(2) - thist(1)); v2z = [v2z, v2z(end)];  % m/s


y0 = [y1z(1); y2z(1)] - thist(1) * [v1z(1); v2z(1)]; % m

% x = [y1z'; y2z'; v1z'; v2z'];

% ψ0 = atan2(vr2,vr1) and vr = ||vr||.

psi0 = atan2(v2z(1), v1z(1)); % rad
psidot0 = 0; % rad/s
vr0 = sqrt(v2z(1)^2 + v1z(1)^2); % m/s

x = [psi0; y0(1); y0(2); psidot0; vr0]; % rad; m; m; rad/s; m/s


za = Rainv' * [rhoahist; rhobhist]; % measurements, m; m



%% plotting
close all

figure
plot(y1z, y2z, 'o')
hold on
plot(la, 0, 'ro')
plot(lb, 0, 'bo')
plot(y0(1), y0(2), 'x')

grid on
axis equal





% fprintf('5: \n\tk_q0: %f\n\tPmd: %f\n', k_q0, Pmd)

%%
% z = ha([y1z; y2z], Rainv, la, lb);

function z = ha(x, Rainv, la, lb)

z = Rainv' * ...
    [vecnorm(x(2:3, :) - [la; 0]); ...
     vecnorm(x(2:3, :) - [lb; 0])]; % m; m

end % function


% function J = 






