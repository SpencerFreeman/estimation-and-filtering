function [xhathist,Phist,sigmahist,enuhist] = ...
    particle_smoother(zkhist, xhat0, P0, Q, R, Ns)

n = length(zkhist); % samples
nx = length(xhat0);
nv = size(Q, 1);

t = 0; % s
xhat = xhat0; % initial state estimate
phat = P0; % initial state covariance
ev = 0;

ts = nan(1, n);
xhats = nan(nx, n);
phats = nan(nx * nx, n);
evs = nan(1, n);

Rinv = inv(R);
Svj = chol(Q)';

for k = 1:n

    ts(k) = t;
    xhats(:, k) = xhat;
    phats(:, k) = phat(:); % unwrap to column vector
    % evs(i) = ev;

    wtil = nan(1, Ns);
    chis = nan(nx, Ns);
    for i = 1:Ns

        chi = chol(P0)'*randn(nx, 1) + xhat0; % initial particle
        for j = 1:k

            vss = Svj * randn(nv, 1);
            chi = f_class_example(j, chi, vss); % propagate
            dz = zkhist(j) - h_class_example(chi);

        end

        wtil(i) = exp( -.5*sum(dz.*Rinv*dz) );
        chis(:, i) = chi; % chi(k)

    end
    w = wtil / sum(wtil);

    xhat = sum(w .* chis, 2); % compute a posteriori state estimate
    phat = zeros(nx);
    for i = 1:Ns
        phat = phat + w(i) * (chis(:, i) - xhat)*(chis(:, i) - xhat)'; % compute a posteriori error covariance matrix
    end % for


end % for

% record the final filter outputs
ts(n) = t;
xhats(:, n) = xhat;
phats(:, n) = phat(:); % unwrap to column vector

xhathist = xhats';
Phist = [];
sigmahist = [];
enuhist = [];


end % function






