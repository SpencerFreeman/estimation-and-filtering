function [J,delxgn,dJdalpha,d2Jdalpha2,P,dJdx] = ...
                     jdxgnmissle(xguess,thist,zhist,lradar,...
                                 sigmarho,sigmatheta,idelx)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function is for the missle tracking problem.  
%  It computes the cost, the Gauss-Newton step, the
%  first and (approximate) second derivatives of the cost with respect
%  to the step size, alpha, the approximate estimation error
%  covariance matrix (assuming that the optimal solution
%  has been reached for the nonlinear least-squares problem), and
%  the first derivative of the cost with respect to x.
%  
%
%  Inputs:
%
%    xguess      = [y10guess;v10guess;y20guess;v20guess], the current
%                guess of the vector of initial conditions of the 
%                missle at time 0.  y10guess is the guessed initial 
%                horizontal position in meters, v10guess is the guessed 
%                initial horizontal velocity in m/sec, y20guess is the 
%                guessed initial vertical position in meters, and 
%                v20guess is the guessed initial vertical velocity
%                in m/sec.
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
%    idelx       A flag that tells whether (idelx = 1) or not (idelx = 1)
%                to compute the other outputs besides J.
%  
%  Outputs:
%
%    J           The nonlinear least-squares cost.  Note that this
%                is divided by 2 compared to the cost given in
%                Bar-Shalom.
%
%    delxgn      The Gauss-Newton perturbation to xguess that is
%                supposed to yield the optimal x.
%
%    dJdalpha    = dJ(xguess+alpha*delxgn)/dalpha evaluated at
%                alpha = 0.
%
%    d2Jdalpha2  = d2J(xguess+alpha*delxgn)/dalpha2 evaluated at
%                alpha = 0, except that the Hessian matrix d2Jdx2 
%                used in this calculation is only the Gauss-Newton 
%                approximation of the cost function Hessian.
%
%    P           The 4x4 Gauss-Newton approximation of the estimation
%                error covariance.  It is the inverse of the approximate
%                Hessian of J.
%
%    dJdx        The 4x1 first partial derivative of J with respect
%                to x.  dJdx(i,1) = dJ/dx(i).
%

%
%  Get the number of radar measurements.
%
   k = size(thist,1);
%
%  Loop through the measurement times and set up the error vector
%  and the large measurement Jacobian matrix.  Re-normalize
%  so as to deal with a nonlinear measurement equation that has
%  a measurement error with the identity matrix for its covariance.
%
   delzabigvec = zeros(2*k,1);
   if idelx == 1
      Habigmat = zeros(2*k,4);
   end
   idumvec = [1:2]';
   Rainv = diag([(1/sigmarho);(1/sigmatheta)]);
   for j = 1:k
      [hjmod,Hjmod] = hjmissle(xguess,thist(j,1),lradar,idelx,0);
      hjmeas = zhist(j,:)';
      delzabigvec(idumvec,1) = Rainv*(hjmeas - hjmod);
      if idelx == 1
         Habigmat(idumvec,:) = Rainv*Hjmod;
      end
      idumvec = idumvec + 2;
   end
%
%  Compute the cost.
%
   J = 0.5*(delzabigvec'*delzabigvec);
%
%  Return if only the cost is needed.
%
   if idelx == 0
      return
   end
%
%  Compute the Gauss-Newton search direction.  Use QR factorization.
%
   [Qb,Rb] = qr(Habigmat,0);
   delzb1bigvec = (Qb')*delzabigvec;
   Rbinv = inv(Rb);
   delxgn = Rbinv*delzb1bigvec;
   P = Rbinv*(Rbinv');
   dJdx = (Habigmat')*delzabigvec;
   dJdalpha = (dJdx')*delxgn;
   dum = Rb*delxgn;
   d2Jdalpha2 = (dum')*dum;