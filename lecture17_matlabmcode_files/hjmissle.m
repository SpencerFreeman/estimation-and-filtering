function [hj,Hj,d2hjdx2] = hjmissle(x,tj,lradar,i1stdrv,i2nddrv)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function gives the measurement function h(x(j),tj) and its 
%  first and second derivatives with respect to x(j), Hj = dhj/dx and
%  d2hjdx2 = dHj/dx. It is for use in the nonlinear least-squares problem 
%  that does missle tracking.  The measurements are range
%  and elevation angle as measured at a radar with a known location.
%
%  This version assumes a flat-Earth gravity field, no
%  aerodynamic drag force, and no other disturbances acting on
%  the missle.
%  
%
%  Inputs:
%
%    x        = [y10;v10;y20;v20], the vector of initial conditions
%             of the missle at time 0.  y10 is the initial horizontal
%             position in meters, v10 is the initial horizontal
%             velocity in m/sec, y20 is the initial vertical position
%             in meters, and v20 is the initial vertical velocity
%             in m/sec.
%
%    tj       The time in seconds of the radar measurement.
%
%    lradar   The horizontal position of the radar station, in meters.
%             The radar station is assumed to be located at zero
%             vertical position.
%
%    i1stdrv  A flag that tells whether (i1stdrv = 1) or not 
%             (i1stdrv = 0) Hj needs to get computed.  Note that
%             second derivatives will not get computed if first
%             derivatives do not get computed.
%
%    i2nddrv  A flag that tells whether (i2nddrv = 1) or not 
%             (i2nddrv = 0) d2hjdx2 needs to get computed.  Note that
%             second derivatives will not get computed if first
%             derivatives do not get computed.  Thus, i2nddrv
%             will be ignored if i1stdrv = 0.
%  
%  Outputs:
%
%    hj       = [rho_j;theta_j], the 2x1 radar output vector.  rho_j
%             is the measured distance from the radar to the missle
%             at time tj given in meters.  theta_j is the measured
%             elevation angle of the missle above the horizon
%             as measured by the radar, in radians.
%
%    Hj       = dhj/dx.  Hj is a 2x4 matrix Hj(1,i) is the derivative
%             of rho_j with respect to x(i).  Hj(2,i) is the derivative
%             of theta_j with respect to x(i).  This output will
%             be needed to to Newton's method or to do the Gauss
%             Newton method.
%
%    d2hjdx2  = d2hj/dx2.  d2hjdx2 is a 2x4x4 array.  d2hjdx2(1,m,n)
%             is the second derivative of rho_j with respect to x(m)
%             and x(n).  d2hjdx2(2,m,n) is the second derivative of 
%             theta_j with respect to x(m) and x(n).  This output will
%             be needed to do Newton's method.
%

%
%  Set up output arrays as needed.
%
   hj = zeros(2,1);
   if i1stdrv == 1
      Hj = zeros(2,4);
      if i2nddrv == 1
         d2hjdx2 = zeros(2,4,4);
      else
         d2hjdx2 = [];
      end
   else
      Hj = [];
      d2hjdx2 = [];
   end
%
%  Compute the position of the missle relative to the radar.
%  Use g = 9.81 m/sec^2 for the acceleration of gravity.
%
   dely1 = lradar - (x(1) + tj*x(2));
   dely2 = x(3) + tj*x(4) - 0.5*9.81*(tj^2);
%
%  Compute the hj outputs.
%
   dely1sq = dely1^2;
   dely2sq = dely2^2;
   hj(1,1) = sqrt(dely1sq + dely2sq);
   hj(2,1) = atan2(dely2,dely1);
%
%  Return if neither first derivatives nor second derivatives
%  need to be calculated.
%
   if i1stdrv == 0
      return
   end
%
%  Calculate the first derivatives.  Use analytic formulas.
%
   ddely1_dx1 = -1;
   ddely1_dx2 = -tj;
   ddely2_dx3 = 1;
   ddely2_dx4 = tj;
   one_over_hj1 = 1/hj(1,1);
   Hj(1,1) = one_over_hj1*dely1*ddely1_dx1;
   Hj(1,2) = one_over_hj1*dely1*ddely1_dx2;
   Hj(1,3) = one_over_hj1*dely2*ddely2_dx3;
   Hj(1,4) = one_over_hj1*dely2*ddely2_dx4;
%
   one_over_hj1sq = one_over_hj1^2;
   Hj(2,1) = -one_over_hj1sq*dely2*ddely1_dx1;
   Hj(2,2) = -one_over_hj1sq*dely2*ddely1_dx2;
   Hj(2,3) =  one_over_hj1sq*dely1*ddely2_dx3;
   Hj(2,4) =  one_over_hj1sq*dely1*ddely2_dx4;
%
%  Return if second derivatives do not need to be calculated.
%
   if i2nddrv == 0
      return
   end
%
%  Calculate the second derivatives.  Use analytic formulas.
%
%  Calculate the 2nd derivatives of hj(1,1).
%
   one_over_hj1cu = one_over_hj1sq*one_over_hj1;
   dely1dely2 = dely1*dely2;
   d2hjdx2(1,1,1) = one_over_hj1cu*dely2sq*(ddely1_dx1^2);
   d2hjdx2(1,1,2) = one_over_hj1cu*dely2sq*(ddely1_dx1*ddely1_dx2);
   d2hjdx2(1,2,1) = d2hjdx2(1,1,2);
   d2hjdx2(1,1,3) = - one_over_hj1cu*dely1dely2*(ddely1_dx1*ddely2_dx3);
   d2hjdx2(1,3,1) = d2hjdx2(1,1,3);
   d2hjdx2(1,1,4) = - one_over_hj1cu*dely1dely2*(ddely1_dx1*ddely2_dx4);
   d2hjdx2(1,4,1) = d2hjdx2(1,1,4);
   d2hjdx2(1,2,2) = one_over_hj1cu*dely2sq*(ddely1_dx2^2);
   d2hjdx2(1,2,3) = - one_over_hj1cu*dely1dely2*(ddely1_dx2*ddely2_dx3);
   d2hjdx2(1,3,2) = d2hjdx2(1,2,3);
   d2hjdx2(1,2,4) = - one_over_hj1cu*dely1dely2*(ddely1_dx2*ddely2_dx4);
   d2hjdx2(1,4,2) = d2hjdx2(1,2,4);
   d2hjdx2(1,3,3) = one_over_hj1cu*dely1sq*(ddely2_dx3^2);
   d2hjdx2(1,3,4) = one_over_hj1cu*dely1sq*(ddely2_dx3*ddely2_dx4);
   d2hjdx2(1,4,3) = d2hjdx2(1,3,4);
   d2hjdx2(1,4,4) = one_over_hj1cu*dely1sq*(ddely2_dx4^2);
%
%  Calculate the 2nd derivatives of hj(2,1).
%
   one_over_hj14th = one_over_hj1cu*one_over_hj1;
   dely1sqmdely2sq = dely1sq - dely2sq;
   d2hjdx2(2,1,1) = 2*one_over_hj14th*dely1dely2*(ddely1_dx1^2);
   d2hjdx2(2,1,2) = 2*one_over_hj14th*dely1dely2*(ddely1_dx1*ddely1_dx2);
   d2hjdx2(2,2,1) = d2hjdx2(2,1,2);
   d2hjdx2(2,1,3) = -one_over_hj14th*dely1sqmdely2sq*...
                                 (ddely1_dx1*ddely2_dx3);
   d2hjdx2(2,3,1) = d2hjdx2(2,1,3);
   d2hjdx2(2,1,4) = -one_over_hj14th*dely1sqmdely2sq*...
                                 (ddely1_dx1*ddely2_dx4);
   d2hjdx2(2,4,1) = d2hjdx2(2,1,4);
   d2hjdx2(2,2,2) = 2*one_over_hj14th*dely1dely2*(ddely1_dx2^2);
   d2hjdx2(2,2,3) = -one_over_hj14th*dely1sqmdely2sq*...
                                 (ddely1_dx2*ddely2_dx3);
   d2hjdx2(2,3,2) = d2hjdx2(2,2,3);
   d2hjdx2(2,2,4) = -one_over_hj14th*dely1sqmdely2sq*...
                                 (ddely1_dx2*ddely2_dx4);
   d2hjdx2(2,4,2) = d2hjdx2(2,2,4);
   d2hjdx2(2,3,3) = -2*one_over_hj14th*dely1dely2*(ddely2_dx3^2);
   d2hjdx2(2,3,4) = -2*one_over_hj14th*dely1dely2*(ddely2_dx3*ddely2_dx4);
   d2hjdx2(2,4,3) = d2hjdx2(2,3,4);
   d2hjdx2(2,4,4) = -2*one_over_hj14th*dely1dely2*(ddely2_dx4^2);