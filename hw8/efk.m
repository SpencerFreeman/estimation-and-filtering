function [xhathist,Phist,sigmahist,enuhist] = ...
    efk(zkhist, xhat0, P0, Q, R)

n = length(zkhist); % samples
nx = length(xhat0);
nv = size(Q, 1);
thist = 1:n;

% EKF
t = 0; % s
xhat = xhat0; % initial state estimate
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
    % tkp1 = thist(i); % s

    % [fprinted, dfprinted_dxk, dfprinted_dvk] = ...
    %     c2dnonlinear(xhat, [], [0; 0], t, tkp1, nRK, fscriptname, true);

    xbar = f_class_example(i, xhat, 0);
    F = 2*sec(xhat)^2; % df / dxk
    GAMMA = 1;

    pbar = F * phat * F' + GAMMA * Q * GAMMA';
    % t = tkp1;

    % measurement update
    zbar = h_class_example(xbar);
    H = 1 + 2*xbar + 3*xbar^2; % dh /dx

    z = zkhist(i);
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

xhathist = xhats';
Phist = reshape(phats, nx, nx, n);
sigmahist = [];
enuhist = [];



end % function

























