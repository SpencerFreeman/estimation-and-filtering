function [sa,dsadx,d2sadx2] = safunct(x)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function computes
%                      _
%                     |  sin(x)/x   if x ~= 0
%            sa(x) = -|
%                     |_ 1          if x == 0
%  
%  and its first 2 derivatives with respect to x.
%
%  Use a series approximation for sin(x) if x is small.
%  
   if abs(x) <= 0.02
      xsq = x^2;
      twox = 2*x;
      d2sadx2 = 0;
      dsadx = 0;
      sa = 1;
      jjfac1 = 9;
      jjfac2 = 8;
      for jj = 1:4
         denfac = 1/(jjfac1*jjfac2);
         d2sadx2 = - (d2sadx2*xsq + 2*dsadx*twox + sa*2)*denfac;
         dsadx = - (dsadx*xsq + sa*twox)*denfac;
         sa = 1 - (sa*xsq)*denfac;
         jjfac1 = jjfac1 - 2;
         jjfac2 = jjfac2 - 2;
      end
   else
      sinx = sin(x);
      cosx = cos(x);
      oneox = 1/x;
      sa = sinx*oneox;
      dsadx = (cosx - sa)*oneox;
      d2sadx2 = (-sinx - 2*dsadx)*oneox;
   end