function [xhathist,Phist,sigmahist,enuhist] = ...
    class_example_ukf(zhist,xhat0,P0, Q, R)


n = size(zhist, 1);
nx = length(xhat0);

nv = size(Q, 1);
nz = size(zhist, 2);


t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
xhats_ukf = nan(nx, n);
phats_ukf = nan(nx * nx, n);
evs = nan(1, n);

Sv = chol(Q)';
alpha = 1;
beta = 2;
K = 1;
% α = 1, κ = 1, and β = 2
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

        fprinted = f_class_example(i, chi(:, j), vss(:, j));
        chi_bar(:, j) = fprinted;
        h = h_class_example(chi_bar(:, j));
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


xhathist = xhats_ukf';
Phist = [];
sigmahist = [];
enuhist = [];

end % function











