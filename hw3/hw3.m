clear;clc;close all

% psi = x(1); % rad
% y10 = x(2); % m
% y20 = x(3); % m
% psid = x(4); % rad/s
% vr  = x(5); % m

idispflag = true; % print outs

load('radarmeasdata_cart.mat')

rhobhist = rhobhist';
rhoahist = rhoahist';
thist = thist';
k = length(thist); % number of samples

nz = 2; % number of measurements
nx = 5; % number of states

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


psi0 = atan2(v2z(1), v1z(1)); % rad
psidot0 = 0; % rad/s
vr0 = sqrt(v2z(1)^2 + v1z(1)^2); % m/s

x0 = [psi0; y0(1); y0(2); psidot0; vr0]; % rad; m; m; rad/s; m/s


z = reshape( Rainv' * [rhoahist; rhobhist], nz*k, 1); % measurements stacked for batch processing, m; m




%% processing loop

xguess = x0;
alpha = 1;

[delx, J, dJdx, dJdalpha, d2Jdalpha2, P] = ...
    newton_gauss(xguess, thist, z, Rainv, la, lb);

delJpred = dJdalpha + .5*d2Jdalpha2; %  Predict the change in cost if a step size of alpha = 1 is taken.

%
%  Prepare some quantities for use in controlling the Gauss-Newton
%  iterations.
%
testdone = 0;
niteration = 0;
iaflag = 0;

tic
while testdone == 0

    alpha = 1;
    xguessnew = xguess + alpha*delx;
    [~, J_new, dJdx, dJdalpha, d2Jdalpha2, P] = ...
        newton_gauss(xguess + alpha * delx, ...
        thist, z, Rainv, la, lb);
    %
    %  Do step size halving if necessary in order to force a decrease
    %  in the cost.
    %
    nalphahalf = 0;
    while J_new >= J
        nalphahalf = nalphahalf + 1;
        if nalphahalf > 50
            iaflag = 1;
            break
        end % if
        alpha = 0.5*alpha;
        xguessnew = xguess + alpha*delx;
        [~, J_new, dJdx, dJdalpha, d2Jdalpha2, P] = ...
            newton_gauss(xguess + alpha * delx, ...
            thist, z, Rainv, la, lb);
    end % while


    % termination conditions
    if iaflag == 1
        itermflag = 1;
        break
    end % if
    xguess = xguessnew;
    Jold = J;
    delJold = J_new - J;
    delJpredold = delJpred;

    [delx, J_new, dJdx, dJdalpha, d2Jdalpha2, P] = ...
        newton_gauss(xguess, ...
        thist, z, Rainv, la, lb);

    delJpred = dJdalpha + .5*d2Jdalpha2;
    delJsizetest = abs(delJpred) < 1.e-13*(1 + J);
    delxsizetest = norm(delx) < 1.e-09*(1 + norm(xguess));
    alphatest = alpha == 1;
    delJratiotest = abs((delJold/delJpredold) - 1) < 0.01;
    if delJsizetest & delxsizetest
        testdone = 1;
    end % if
    if alphatest & delJratiotest & delxsizetest
        testdone = 1;
    end % if
    niteration = niteration + 1;
    if testdone == 0
        if niteration >= 100
            itermflag = 2;
            testdone = 1;
        end % if
    end % if


    if idispflag == 1
        disp([' At iteration ',int2str(niteration),' alpha = ',...
            num2str(alpha),', Jnew = ',num2str(J),', Jold = ',...
            num2str(Jold),', and norm(delxnew) = ',...
            num2str(norm(delx)),'.'])
    end % if




end % for
toc

xest = xguess;
Jopt = J;












%% plotting
close all

figure
plot(y1z, y2z, 'o')
hold on
plot(la, 0, 'ro')
plot(lb, 0, 'bo')
plot(y0(1), y0(2), 'x')
plot(xest(2), xest(3), 'oc')

grid on
axis equal





fprintf('5: \n\tXest: %f %f %f %f %f\n\tJopt: %f\n', xest, Jopt)

%%
function [delx, J, dJdx, dJdalpha, d2Jdalpha2, P] = ...
    newton_gauss(xguess, thist, z, Rainv, la, lb)


    H = dhdx_a(xguess, thist, Rainv, la, lb);
    hg = h_a(xguess, thist, Rainv, la, lb);

    [Qb, Rb] = qr(H, 0);

    Rbinv = inv(Rb);

    innov = z - hg; % estimator discrepancy

    delx = Rbinv*Qb'*innov; % x step direction

    J = 0.5*(innov'*innov); % cost function

    P = Rbinv*Rbinv';
    dJdx = H'*innov;
    dJdalpha = dJdx'*delx;
    dum = Rb*delx;
    d2Jdalpha2 = dum'*dum;

end % function




function za = h_a(x, t, Rainv, la, lb)

nz = 2; % number of measurements
k = length(t); % numnber of samples

sa = nan(1, k);
for i = 1:k
    [sa(i), ~, ~] = safunct(.5*x(4)*t(i));
end % for

ys = [x(2) + x(5)*t.*sa.*cos(x(1) + .5*x(4)*t); ... % m
      x(3) + x(5)*t.*sa.*sin(x(1) + .5*x(4)*t)]; % m

za = reshape( Rainv' * ...
    [vecnorm( ys - [la; 0] ); ...
     vecnorm( ys - [lb; 0] )], nz*k, 1);

end % function


function Ha = dhdx_a(x, t, Rainv, la, lb)

k = length(t);
nx = length(x);
nz = 2;

Ha = nan(k*nz, nx);

for i = 1:k

    Ha((i - 1)*nz + (1:nz), :) = ...
        Rainv' * ...
        [drhodx(x, t(i), la); ...
         drhodx(x, t(i), lb)];

end % for


end % function


function dh = drhodx(x, t, l)

theta = x(1) + .5*x(4)*t; % rad
x5t = x(5)*t; % m
ct = cos(theta);
st = sin(theta);

[sa, dsadx, ~] = safunct(.5*x(4)*t);
dsadx4 = .5*t*dsadx;

k0 = .5*( (x(2) + x5t*sa*ct - l)^2 + ...
    (x(3) + x5t*sa*st    )^2 )^(-1/2);

k1 = 2*(x(2) + x5t*sa*ct - l);
k2 = 2*(x(3) + x5t*sa*st - l);

dhdx1 = ...
    ( x5t*sa*(-1)*st * k1 + ...
    x5t*sa*     ct * k2) * k0;

dhdx2 = (1)*k1*k0;

dhdx3 = (1)*k2*k0;

dhdx4 = ...
    ( x5t*( dsadx4*ct + sa*.5*t*(-1)*st ) * k1 + ...
    x5t*( dsadx4*st + sa*.5*t*     ct ) * k2 ) * k0;

dhdx5 = ( t*sa*ct * k1 + t*sa*st * k2 ) * k0;

dh = [dhdx1, dhdx2, dhdx3, dhdx4, dhdx5];

end % function






















