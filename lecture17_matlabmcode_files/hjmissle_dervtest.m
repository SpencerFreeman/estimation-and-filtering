%hjmissle_dervtest.m
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This matlab script tests the derivative calculations in hjmissle.m
%  by performing finite-difference calculations of the derivatives
%  and comparing them to the analytic derivatives.
%
   x = [101.e+03;3.0e+03;205.4;2.5e+03];
   tj = 230.9;
   lradar = 1.2e+06;
   [hj,Hj,d2hjdx2] = hjmissle(x,tj,lradar,1,1);
   Hj_num = zeros(2,4);
   d2hjdx2_num = zeros(2,4,4);
%
%  Appropriate finite-difference intervals for the numerical
%  differentiation are problem dependent.  One often must "play"
%  with these in order to get reasonable values. Values that are
%  too large will yield problems with the finite-difference
%  formula not being a good approximation because of higher-order
%  terms.  Values that are too small with yield inaccuracies due to
%  computer round-off errors.
%  
   delxvec = [10;.01;10;.01];
%
%  Use central differencing to get numerical approximations
%  to the derivatives.
%
   for jj = 1:4
      xdum = x;
      xdum(jj,1) = x(jj,1) + delxvec(jj,1);
      [hjp,Hjp] = hjmissle(xdum,tj,lradar,1,0);
      xdum(jj,1) = x(jj,1) - delxvec(jj,1);
      [hjm,Hjm] = hjmissle(xdum,tj,lradar,1,0);
      delxfac = 1/(2*delxvec(jj,1));
      Hj_num(:,jj) = (hjp - hjm)*delxfac;
      d2hjdx2_num(:,:,jj) = (Hjp - Hjm)*delxfac;
   end
   clear jj xdum hjp Hjp hjm Hjm
%
%  Print out checks of the errors.  
%
%  These are conservative.  These should be small, on the order
%  of 1.e-05 or smaller, otherwise there almost definitely is a 
%  problem.
%
   normHjtest = norm(Hj - Hj_num)/norm(Hj_num)
   dum = d2hjdx2 - d2hjdx2_num;
   dum = dum(:);
   normd2hjdx2test = norm(dum)/norm(d2hjdx2_num(:))
%
%  These are less conservative.  These ought to be small, on the order
%  of 1.e-05 or smaller, but if they are not small, then look closer
%  at the arrays in order to see why these test outputs
%  are not small before deciding if there is a problem.  
%  The eps^3 term gets added to avoid division by zero when
%  the correct term is supposed to be zero.  Except for the eps^3
%  terms this tests the relative error in each element of the
%  array and gives the worst-case relative error.
%
   dum = Hj - Hj_num;
   dum = dum(:);
   dum1 = Hj_num(:);
   maxHjtest = max(abs(dum)./(abs(dum1) + eps^3))
   dum = d2hjdx2 - d2hjdx2_num;
   dum = dum(:);
   dum1 = d2hjdx2_num(:);
   maxd2hjdx2test = max(abs(dum)./(abs(dum1) + eps^3))
   clear dum dum1
%
%  If there is a problem with either test, then
%  the first thing to try is to see if an increase or a decrease of
%  delxvec (by a factor of 10) has a positive impact on the sizes
%  of all 4 test quantities.
%  