function [J,delxnewton,dJdalpha,d2Jdalpha2,P,dJdx,d2Jdx2,iflagHess] = ...
                     jdxnewtonmissle(xguess,thist,zhist,lradar,...
                                     sigmarho,sigmatheta,idelx)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function is for the missle tracking problem.  
%  It computes the cost, the Newton step, the
%  first and second derivatives of the cost with respect
%  to the step size, alpha, the approximate estimation error
%  covariance matrix (assuming that the optimal solution
%  has been reached for the nonlinear least-squares problem), and
%  the first and second derivatives of the cost with respect to x.
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
%    delxnewton  The Newton perturbation to xguess that is
%                supposed to yield the optimal x.
%
%    dJdalpha    = dJ(xguess+alpha*delxgn)/dalpha evaluated at
%                alpha = 0.
%
%    d2Jdalpha2  = d2J(xguess+alpha*delxgn)/dalpha2 evaluated at
%                alpha = 0.
%
%    P           The 4x4 Newton approximation of the estimation
%                error covariance.  It is the inverse of the approximate
%                Hessian of J.  This will be the inverse of the
%                exact d2Jdx2 Hessian matrix if it is positive 
%                definite.
%
%    dJdx        The 4x1 first partial derivative of J with respect
%                to x.  dJdx(i,1) = dJ/dx(i).
%
%    d2Jdx2      The 4x4 second partial derivative of J with respect
%                to x.  d2Jdx2(m,n) = d2J/dx(m)dx(n).
%
%    iflagHess   A flag that tells whether (iflagHess = 1) or not
%                (iflagHess = 0) the Hessian of the cost function
%                had to get modified by adding a positive
%                definite matrix to it in order to get it to be
%                positive definite so that the Newton calculation
%                would produce a direction in which a cost decrease
%                could be guaranteed to occur.

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
   V2ndO = zeros(4,4);
   for j = 1:k
      [hjmod,Hjmod,d2hjdx2] = hjmissle(xguess,thist(j,1),lradar,idelx,1);
      hjmeas = zhist(j,:)';
      delzabigvec(idumvec,1) = Rainv*(hjmeas - hjmod);
      if idelx == 1
         Habigmat(idumvec,:) = Rainv*Hjmod;
         dum = d2hjdx2(1,:,:);
         dum1 = zeros(4,4);
         dum1(:) = dum(:);
         dum1 = dum1*(Rainv(1,1)*delzabigvec(1,1));
         V2ndO = V2ndO - dum1;
         dum = d2hjdx2(2,:,:);
         dum1(:) = dum(:);
         dum1 = dum1*(Rainv(2,2)*delzabigvec(2,1));
         V2ndO = V2ndO - dum1;
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
%  Compute the Newton search direction.  Use QR factorization.
%
   [Qb,Rb] = qr(Habigmat,0);
   delzb1bigvec = (Qb')*delzabigvec;
   Rbinv = inv(Rb);
%
%  The following calculations incorporate the 2nd derivative
%  terms of h(x) into the Hessian.  They do this in a way that
%  preserves the use of square-root calculations.  Also, if 
%  the resulting cost Hessian is not positive definite, 
%  then a positive definite increment gets added for
%  purposes of computing the Newton search direction so that
%  it can be guaranteed to be a descent direction of the cost
%  function.
%
   Dmat = eye(4) + (Rbinv')*(V2ndO*Rbinv);
%
%  Compute the Newton search direction.  Use QR factorization.
%
   evals = eig(Dmat);
   elim = 1.e-12*max(abs(evals));
   emin = min(evals);
   iflagHess = 0;
   if emin < elim
      Dmat = Dmat + eye(4)*(elim - emin);
      iflagHess = 1;
   end
   Rmat = chol(Dmat);
   Rmat_inv = inv(Rmat);
   Rtotinv = Rbinv*Rmat_inv;
   delxnewton = Rtotinv*((Rmat_inv')*delzb1bigvec);
   P = Rtotinv*(Rtotinv');
   dJdx = (Habigmat')*delzabigvec;
   d2Jdx2 = (Rb')*Rb + V2ndO;
   dJdalpha = (dJdx')*delxnewton;
   d2Jdalpha2 = (delxnewton')*(d2Jdx2*delxnewton);