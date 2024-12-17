clear;clc;close all

do_ukf = true;

cart_EKF_meas


b = 0.1; % m
t_steer = 0.25; % sec
t_speed = 0.60; % sec
xbar_5 = 2.1; % m/sec
la = -1; % m
lb = 1; % m

% The radar measurement noise standard deviations are
sig_rho_a = 0.002; % m
sig_rho_b = 0.002; % m

% The continuous-time process noise intensities are
q_til_steer = 0.10; % rad^2/sec
q_til_speed = 21.25; % m^2/sec^3

%% develop initial guess of x

rhobhist = zhist(:, 1); % m
rhoahist = zhist(:, 2); % m

dt = thist(2) - thist(1); % s

y1z = (rhobhist.^2 - rhoahist.^2 - lb^2 + la^2) / 2 / (la - lb); % m
y2z = sqrt(rhoahist.^2 - (y1z - la).^2); % m

v1z = diff(y1z) / dt; v1z = [v1z; v1z(end)]; % m/s
v2z = diff(y2z) / dt; v2z = [v2z; v2z(end)];  % m/s

y0 = [y1z(1); y2z(1)]; % - thist(1) * [v1z(1); v2z(1)]; % m

psi0 = atan2(v2z(1), v1z(1)); % rad
psidot0 = 0; % rad/s
vr0 = sqrt(v2z(1)^2 + v1z(1)^2); % m/s

x0 = [psi0; y0(1); y0(2); psidot0; vr0]; % initial guess, rad; m; m; rad/s; m/s
% x0(2:3) = [-1;2];
% x0 = [1, .4, 1.5, 0, 6.7]';

P0 = diag([ ...
    2*sig_rho_a*sig_rho_b/vr0^2/dt^2, ...
    sig_rho_a*sig_rho_b, ...
    sig_rho_a*sig_rho_b, ...
    q_til_steer*t_steer/2, ...
    2*sig_rho_a*sig_rho_b/dt^2]);

%%

fscriptname = 'fscript_cart';

Q = diag([q_til_steer, q_til_speed]) / dt;
R = diag([sig_rho_a^2, sig_rho_b^2]); % m^2

nRK = 60; % RK steps

n = length(thist) + 1;
nx = length(x0);
nv = 2;
nz = 2;

%% UKF
t = 0; % s
xhat = x0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
xhats_ukf = nan(nx, n);
phats_ukf = nan(nx * nx, n);
evs = nan(1, n);

if do_ukf
    Sv = chol(Q)';
    alpha = 1;
    beta = 2;
    K = 3 - nx - nv;
    lam = alpha^2 * (nx + nv + K) - (nx + nv);
    sqrt_nx_nv_lam = sqrt(nx + nv + lam);
    W_m_0 = lam/(nx + nv + lam);
    W_m_i = 1/2/(nx + nv + lam);
    W_c_0 = lam/(nx + nv + lam) + 1 - alpha^2 + beta;
    W_c_i = 1/2/(nx + nv + lam);

    for i = 1:(n - 1)

        ts(i) = t;
        xhats_ukf(:, i) = xhat;
        phats_ukf(:, i) = phat(:); % unwrap to column vector
        % evs(i) = ev;

        % propagate
        tkp1 = thist(i); % s

        Sx = chol(phat)';

        chi = nan(nx, 2*nx + 2*nv);
        vss = nan(nv, 2*nx + 2*nv);
        for i_sig = 1:nx
            chi(:, i_sig) = xhat + sqrt_nx_nv_lam * Sx(:, i_sig);
            vss(:, i_sig) = zeros(nv, 1);
        end
        for i_sig = (nx + 1):(2*nx)
            chi(:, i_sig) = xhat - sqrt_nx_nv_lam * Sx(:, i_sig - nx);
            vss(:, i_sig) = zeros(nv, 1);
        end
        for i_sig = (2*nx + 1):(2*nx + nv)
            chi(:, i_sig) = xhat;
            vss(:, i_sig) = sqrt_nx_nv_lam * Sv(:, i_sig - 2*nx);
        end
        for i_sig = (2*nx + nv + 1):(2*nx + 2*nv)
            chi(:, i_sig) = xhat;
            vss(:, i_sig) = -sqrt_nx_nv_lam * Sv(:, i_sig - 2*nx - nv);
        end
        chi = [xhat, chi];
        vss = [zeros(nv, 1), vss];

        chi_bar = nan(nx, size(chi, 2));
        zss_bar = nan(nz, size(chi, 2));
        for j = 1:size(chi, 2)

            [fprinted, ~, ~] = ...
                c2dnonlinear(chi(:, j), [], vss(:, j), t, tkp1, nRK, fscriptname, false);

            chi_bar(:, j) = fprinted;

            [h, ~] = h_cart(chi_bar(:, j), false);
            zss_bar(:, j) = h;

        end % for

        xbar = W_m_0 * chi_bar(:, 1) + ...
            sum(W_m_i * chi_bar(:, 2:end), 2);
        zbar = W_m_0 * zss_bar(:, 1) + ...
            sum(W_m_i * zss_bar(:, 2:end), 2);

        pbar = W_c_0 * (chi_bar(:, 1) - xbar) * (chi_bar(:, 1) - xbar)';
        pxzbar = W_c_0 * (chi_bar(:, 1) - xbar) * (zss_bar(:, 1) - zbar)';
        pzzbar = W_c_0 * (zss_bar(:, 1) - zbar) * (zss_bar(:, 1) - zbar)' + R;

        for j = 2:size(chi_bar, 2)

            pbar =  pbar + W_c_i * (chi_bar(:, j) - xbar) * (chi_bar(:, j) - xbar)';
            pxzbar = pxzbar + W_c_i * (chi_bar(:, j) - xbar) * (zss_bar(:, j) - zbar)';
            pzzbar = pzzbar + W_c_i * (zss_bar(:, j) - zbar) * (zss_bar(:, j) - zbar)' + R;

        end % for

        t = tkp1;

        % measurement update
        z = zhist(i, :)';
        v = z - zbar; % innovation
        pzzbarinv = inv(pzzbar);
        pxzbar_pzzbarinv = pxzbar * pzzbarinv;
        xhat = xbar + pxzbar_pzzbarinv * v;
        phat = pbar - pxzbar_pzzbarinv * pxzbar';

        ev = nan;%v' * Sinv * v; % estimation error statistic

    end % for

    % record the final filter outputs
    ts(n) = t;
    xhats_ukf(:, n) = xhat;
    phats_ukf(:, n) = phat(:); % unwrap to column vector

end % if


%% EKF
t = 0; % s
xhat = x0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
vs = nan(1, n);
xhats = nan(nx, n);
phats = nan(nx * nx, n);
evs = nan(1, n);

for i = 1:(n - 1)

    ts(i) = t;
    xhats(:, i) = xhat;
    phats(:, i) = phat(:); % unwrap to column vector
    evs(i) = ev;

    % propagate
    tkp1 = thist(i); % s

    [fprinted, dfprinted_dxk, dfprinted_dvk] = ...
        c2dnonlinear(xhat, [], [0; 0], t, tkp1, nRK, fscriptname, true);

    xbar = fprinted;
    F = dfprinted_dxk;
    GAMMA = dfprinted_dvk;

    pbar = F * phat * F' + GAMMA * Q * GAMMA';
    t = tkp1;

    % measurement update
    [h, dh_dx] = h_cart(xbar, true);
    zbar = h;
    H = dh_dx;

    z = zhist(i, :)';
    v = z - zbar; % innovation
    S = H * pbar * H' + R; Sinv = inv(S);
    W = pbar * H' * Sinv;

    xhat = xbar + W * v;
    phat = pbar - W * S * W';

    ev = v' * Sinv * v; % estimation error statistic

end % for

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector
evs(n) = ev;





%% plotting
close all

% ground track
f = figure;
f.WindowStyle = 'Docked';
plot(xhats(2, :), xhats(3, :), '*')
hold on
plot(xhats_ukf(2, :), xhats_ukf(3, :), 'o')
plot(x0(2), x0(3), 'rx')

title('Estimated Position of Midpoint Between Rear Wheels of Cart')
legend('EKF', 'UKF', 'Initial Guess')
grid on
axis equal
xlabel('East (m)')
ylabel('North (m)')

% time histories
names = ["Heading (rad)", "East (m)", "North (m)", "Steer Angle (rad)", "Speed (m/s)"];

f = figure;
f.WindowStyle = 'Docked';
for i = 1:nx
    subplot(nx, 1, i)
    plot(ts, xhats(i, :), '*'); hold on; grid on
    plot(ts, xhats_ukf(i, :), 'o')
    ylabel(names(i))
    if i == 1
        title('Estimated State Time Histories')
        legend('EKF', 'UKF')
    end % if
end % for
xlabel('Time (s)')

grid on
























