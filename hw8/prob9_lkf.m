function [xhathist,Phist,sigmahist,enuhist] = ...
    prob9_lkf(zkhist, xhat0, P0, Q, R, Gamma, F, H)

n = length(zkhist);
nx = length(xhat0);
 
t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance

ts = nan(1, n);
xhats = nan(nx, n);
phats = nan(nx * nx, n);
for i = 1:(n - 1)

    ts(i) = t;
    xhats(:, i) = xhat;
    phats(:, i) = phat(:); % unwrap to column vector

    % t = thist(i); % s
    xbar = F * xhat; % propagate state estimate
    pbar = F * phat * F' + Gamma * Q * Gamma'; % propagate state covariance

    zbar = H * xbar; % expected measurement 
    z = zkhist(i); % actual measurement
    v = z - zbar; % filter innovation

    S = H * pbar * H' + R; % expected measurement covariance
    W = pbar * H' * inv(S); % filter gain

    xhat = xbar + W * v; % updated state estimate
    phat = pbar - W * S * W'; % updates state covariance

end

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector


xhathist = xhats';
Phist = [];
sigmahist = [];
enuhist = [];











end % function