%% Solve Nonlinear Least-Squares Problem using the Gauss-Newton Method
% Spencer Freeman, 10/18/2024
% AOE 5784, Estimation and Filtering
%
% This script and functions solves number 5 of problem set 3.
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW3-P5')

%% initialize

idispflag = true; % print outs

load('radarmeasdata_cart.mat') % load measurement data

rhobhist = rhobhist'; % radar b history
rhoahist = rhoahist'; % radar a history
thist = thist';
k = length(thist); % number of samples

nz = 2; % number of measurements
nx = 5; % number of states

la = -1; % radar a y-positon, m
lb = 1; % radar b y-positon, m

sig_rhoa = .005; % radar a STD, m
sig_rhob = .005; % radar b STD, m

R = diag([sig_rhoa^2, sig_rhob^2]); % measurement covariance, m^2
Ra = chol(R); % Cholesky-factor square root, m^2
Rainv = inv(Ra); % 

z = reshape( Rainv' * [rhoahist; rhobhist], nz*k, 1); % stack measurements for batch processing, m

%% develop initial guess of x

y1z = (rhobhist.^2 - rhoahist.^2 - lb^2 + la^2) / 2 / (la - lb); % m
y2z = sqrt(rhoahist.^2 - (y1z - la).^2); % m

v1z = diff(y1z) / (thist(2) - thist(1)); v1z = [v1z, v1z(end)]; % m/s
v2z = diff(y2z) / (thist(2) - thist(1)); v2z = [v2z, v2z(end)];  % m/s

y0 = [y1z(1); y2z(1)] - thist(1) * [v1z(1); v2z(1)]; % m

psi0 = atan2(v2z(1), v1z(1)); % rad
psidot0 = 0; % rad/s
vr0 = sqrt(v2z(1)^2 + v1z(1)^2); % m/s

x0 = [psi0; y0(1); y0(2); psidot0; vr0]; % initial guess, rad; m; m; rad/s; m/s
% x0(2:3) = [-1;2];

%% processing loop

xguess = x0;
alpha = 1;

[delx, J, dJdx, dJdalpha, d2Jdalpha2, P] = ...
    newton_gauss(xguess, thist, z, Rainv, la, lb); % cost at initial guess

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
        newton_gauss(xguessnew, thist, z, Rainv, la, lb);
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
            newton_gauss(xguessnew, thist, z, Rainv, la, lb);
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
        newton_gauss(xguess, thist, z, Rainv, la, lb);

    delJpred = dJdalpha + .5*d2Jdalpha2;
    delJsizetest = abs(delJpred) < 1.e-13*(1 + J);
    delxsizetest = norm(delx) < 1.e-09*(1 + norm(xguess));
    alphatest = alpha == 1;
    delJratiotest = abs((delJold/delJpredold) - 1) < 0.01;
    if delJsizetest && delxsizetest
        testdone = 1;
    end % if
    if alphatest & delJratiotest && delxsizetest
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

h = figure;
h.WindowStyle = 'Docked';
plot(y1z, y2z, 'o')
hold on
plot(la, 0, 'o', 'Color', "#EDB120")
plot(lb, 0, 'o', 'Color', "#EDB120")
plot(x0(2), x0(3), 'x')
plot(xest(2), xest(3), '*', 'Color', "#D95319")

[psihist, y1hist, y2hist] = psiy1y2cart(xest, (0:.01:thist(end))');
plot(y1hist, y2hist, '-', 'Color', "#D95319")

grid on
axis equal
legend('Measurements', 'Radar A', 'Radar B', 'Guess', 'Estimate', '', ...
    'Location', 'SouthEast')

fprintf('\nXguess: %f %f %f %f %f\nXest: %f %f %f %f %f\nJopt: %f\n', x0, xest, Jopt)

%%
function [delx, J, dJdx, dJdalpha, d2Jdalpha2, P] = ...
    newton_gauss(xguess, thist, z, Rainv, la, lb)

[hg, H] = predicted_measurements(xguess, thist, Rainv, la, lb);

[Qb, Rb] = qr(H, 0); % QR factorization
Rbinv = inv(Rb);

innov = z - hg; % estimator-measurement discrepancy
delx = Rbinv*Qb'*innov; % x step direction

J = 0.5*(innov'*innov); % cost function

P = Rbinv*Rbinv';
dJdx = H'*innov;
dJdalpha = dJdx'*delx;
dum = Rb*delx;
d2Jdalpha2 = dum'*dum;

end % function

function [za, Ha] = predicted_measurements(x, t, Rainv, la, lb)

nx = length(x); % number of states
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






















