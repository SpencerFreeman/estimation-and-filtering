function [xest,Jopt,P,itermflag] = ...
    gnestmissle(thist,zhist,lradar,...
    sigmarho,sigmatheta,idispflag)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.
%
%  This function does Gauss-Newton estimation of the initial
%  state of a missle in a missle tracking problem.
%
%
%  Inputs:
%
%    thist       The kx1 vector of radar measurement times.
%
%    zhist       = [rho_1,theta_1;rho_2,theta_2;rho_3,theta_3;...;...
%                rho_k,theta_k], the kx2 time history of the
%                radar measurements.
%
%    lradar      The horizontal position of the radar station, in meters.
%                The radar station is assumed to be located at zero
%                vertical position.
%
%    sigmarho    The measurement error standard deviation of the
%                radar range measurement, in meters.
%
%    sigmatheta  The measurement error standard deviation of the
%                radar elevation measurement, in radians.
%
%    idispflag   A flag that tells whether (idispflag = 1) or not
%                (idispflag = 0) to display interm optimization
%                results during the Gauss-Newton iterations.
%
%  Outputs:
%
%    xest        = [y10est;v10est;y20est;v20est], the final
%                estimate of the vector of initial conditions of the
%                missle at time 0.  y10est is the estimated initial
%                horizontal position in meters, v10est is the estimated
%                initial horizontal velocity in m/sec, y20est is the
%                estimated initial vertical position in meters, and
%                v20est is the estimated initial vertical velocity
%                in m/sec.
%
%    Jopt        The final optimal value of the nonlinear least-squares
%                cost.  Note that this is divided by 2 compared to the
%                cost given in Bar-Shalom.
%
%    P           The 4x4 Gauss-Newton approximation of the estimation
%                error covariance.  It is the inverse of the approximate
%                Hessian of J.
%
%    itermflag   A termination flag that tells whether the
%                procedure terminated well or not.  Its output values
%                have the following interpretations:
%
%                   0   Normal termination at the optimum
%
%                   1   Terminated because more than 50 step
%                       size halvings were required.  This
%                       may be an optimum.
%
%                   2   Terminated because more than 100
%                       Gauss-Newton iterations were required.
%                       This may not be the optimum.
%

%
%  Get the number of radar measurements.
%
k = size(thist,1);
%
%  Make a first guess based on the first and last measurements
%  and a belief that there is zero measurement error.
%
rho_1 = zhist(1,1);
rho_k = zhist(k,1);
theta_1 = zhist(1,2);
theta_k = zhist(k,2);
t_1 = thist(1,1);
t_k = thist(k,1);
y1_1 = lradar - rho_1*cos(theta_1);
y1_k = lradar - rho_k*cos(theta_k);
y2_1 = rho_1*sin(theta_1);
y2_k = rho_k*sin(theta_k);
Hdum = [1 t_1 0 0;0 0 1 t_1;1 t_k 0 0;0 0 1 t_k];
delz = [y1_1;(y2_1 + 0.5*9.81*(t_1^2));y1_k;(y2_k + 0.5*9.81*(t_k^2))];
xguess = Hdum\delz;
%
%  Calculate the initial cost and the initial Gauss-Newton search
%  direction.
%
[J,delxgn,dJdalpha,d2Jdalpha2,P,dJdx] = ...
    jdxgnmissle(xguess,thist,zhist,lradar,...
    sigmarho,sigmatheta,1);
%
%  Predict the change in cost if a step size of alpha = 1 is taken.
%
delJpred = dJdalpha + .5*d2Jdalpha2;
%
%  Decide whether to terminate now.
%
delJsizetest = abs(delJpred) < 1.e-13*(1 + J);
delxsizetest = norm(delxgn) < 1.e-09*(1 + norm(xguess));
itermflag = 0;
if delJsizetest & delxsizetest
    xest = xguess;
    Jopt = J;
    return
end % if
%
%  Prepare some quantities for use in controlling the Gauss-Newton
%  iterations.
%
testdone = 0;
niteration = 0;
iaflag = 0;
%
%  Do one Gauss-Newton iteration per iteration of the following while
%  loop.
%
while testdone == 0
    alpha = 1;
    xguessnew = xguess + alpha*delxgn;
    Jnew = jdxgnmissle(xguessnew,thist,zhist,lradar,...
        sigmarho,sigmatheta,0);
    %
    %  Do step size halving if necessary in order to force a decrease
    %  in the cost.
    %
    nalphahalf = 0;
    while Jnew >= J
        nalphahalf = nalphahalf + 1;
        if nalphahalf > 50
            iaflag = 1;
            break
        end % if
        alpha = 0.5*alpha;
        xguessnew = xguess + alpha*delxgn;
        Jnew = jdxgnmissle(xguessnew,thist,zhist,lradar,...
            sigmarho,sigmatheta,0);
    end % while

    if iaflag == 1
        itermflag = 1;
        break
    end % if
    xguess = xguessnew;
    Jold = J;
    delJold = Jnew - J;
    delJpredold = delJpred;
    [J,delxgn,dJdalpha,d2Jdalpha2,P,dJdx] = ...
        jdxgnmissle(xguess,thist,zhist,lradar,...
        sigmarho,sigmatheta,1);

    delJpred = dJdalpha + .5*d2Jdalpha2;
    delJsizetest = abs(delJpred) < 1.e-13*(1 + J);
    delxsizetest = norm(delxgn) < 1.e-09*(1 + norm(xguess));
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
            num2str(norm(delxgn)),'.'])
    end % if
    
end % while
xest = xguess;
Jopt = J;

















