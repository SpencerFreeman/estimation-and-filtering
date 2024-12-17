function [xhathist,Phist,sigmahist,enuhist] = ...
    particle_filter(zkhist, xhat0, P0, Q, R)

n = length(zkhist); % samples
nx = length(xhat0);
nv = size(Q, 1);
thist = 1:n;

t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
xhats_pft = nan(nx, n);
phats_pft = nan(nx * nx, n);
evs = nan(1, n);

Rinv = inv(R);

Ns = 400; % # of particles

w0 = 1/Ns * ones(1, Ns); % initial weights
chi0 = chol(P0)'*randn(nx, Ns) + xhat0; % initial particles

w = w0;
chi = chi0;
for i = 1:n

    ts(i) = t;
    xhats_pft(:, i) = xhat;
    phats_pft(:, i) = phat(:); % unwrap to column vector
    % evs(i) = ev;

    vss = chol(P0)'*randn(nv, Ns);

    z = zkhist(i);

    log_wtil = nan(1, Ns);
    for j = 1:Ns
        chi(:, j) = f_class_example(i, chi(:, j), vss(:, j)); % propagate to k+1

        log_wtil(j) = log(w(j)) - .5 * (z - h_class_example(chi(:, j)))' * Rinv * (z - h_class_example(chi(:, j)));
        % wtil(j) = w(j) * exp( -.5 * (z - h(chi(:, j)))' * Rinv * (z - h(chi(:, j))) );
    end % for
    log_wtil_max = max(log_wtil);
    wtiltil = exp(log_wtil - log_wtil_max);

    w = wtiltil / sum(wtiltil); % normalized weights

    xhat = sum(w .* chi, 2); % compute a posteriori state estimate
    phat = zeros(nx);
    for j = 1:Ns
        phat = phat + w(j) * (chi(:, j) - xhat)*(chi(:, j) - xhat)'; % compute a posteriori error covariance matrix
    end % for

    % resampling
    c = nan(1, Ns + 1);
    c(1) = 0;
    c(end) = 1 + 10^-10;
    for j = 2:Ns
        c(j) = sum(w(1:j - 1));
    end % for

    chi_new = nan(nx, Ns);
    for l = 1:Ns
        nl = rand;
        ind = find(nl >= c, 1, 'last');
        chi_new(:, l) = chi(:, ind);
    end % for

    chi = chi_new;
    w = w0;

end % for

% record the final filter outputs
ts(n) = t;
xhats_pft(:, n) = xhat;
phats_pft(:, n) = phat(:); % unwrap to column vector

xhathist = xhats_pft';
Phist = [];
sigmahist = [];
enuhist = [];


end % function






